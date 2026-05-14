import { useState, useRef, useEffect, useMemo } from "react";
import ReactMarkdown from "react-markdown";
import remarkGfm from "remark-gfm";
import rehypeHighlight from "rehype-highlight";
import { useLocale } from "../../hooks/useLocale";
import { useStore } from "../../store";

interface TokenUsage {
  prompt_tokens: number;
  completion_tokens: number;
  total_tokens: number;
}

interface ToolCallSeg {
  type: "tool_call";
  tool_name: string;
  args: string;
  result: string;
}

interface TextSeg {
  type: "text";
  content: string;
}

type Segment = TextSeg | ToolCallSeg;

interface Message {
  id: string;
  role: string;
  segments: Segment[];
  isError?: boolean;
  tokenUsage?: TokenUsage;
}

function buildMessagesFromApi(rawMessages: any[]): Message[] {
  const result: Message[] = [];
  let i = 0;
  while (i < rawMessages.length) {
    const m = rawMessages[i];

    if (m.role === "tool") {
      i++;
      continue;
    }

    const segments: Segment[] = [];

    if (m.content) {
      segments.push({ type: "text", content: m.content });
    }

    if (m.tool_calls && Array.isArray(m.tool_calls)) {
      for (const tc of m.tool_calls) {
        let resultText = "";
        let j = i + 1;
        while (j < rawMessages.length && rawMessages[j].role === "tool") {
          if (rawMessages[j].tool_call_id === tc.id) {
            resultText = rawMessages[j].content || "";
            break;
          }
          j++;
        }
        const fn = tc.function || {};
        segments.push({
          type: "tool_call",
          tool_name: fn.name || "",
          args: fn.arguments || "",
          result: resultText,
        });
      }
    }

    let tokenUsage: TokenUsage | undefined;
    const metaUsage = m.metadata?.token_usage;
    if (metaUsage) {
      tokenUsage = {
        prompt_tokens: metaUsage.prompt_tokens || 0,
        completion_tokens: metaUsage.completion_tokens || 0,
        total_tokens: metaUsage.total_tokens || 0,
      };
    }

    result.push({
      id: m.id,
      role: m.role,
      segments,
      tokenUsage,
    });
    i++;
  }
  return result;
}

export function ChatPanel() {
  const { t } = useLocale();
  const { currentSessionId, setSessions, setCurrentSessionId } = useStore();
  const [messages, setMessages] = useState<Message[]>([]);
  const [input, setInput] = useState("");
  const [streaming, setStreaming] = useState(false);
  const [streamingSegments, setStreamingSegments] = useState<Segment[]>([]);
  const messagesEnd = useRef<HTMLDivElement>(null);
  const abortRef = useRef<AbortController | null>(null);
  const skipHistory = useRef(false);
  const pendingUsageRef = useRef<TokenUsage | undefined>(undefined);

  const estimateInputTokens = (text: string) => {
    let cjk = 0;
    for (let i = 0; i < text.length; i++) {
      const c = text.charCodeAt(i);
      if ((c >= 0x4e00 && c <= 0x9fff) || (c >= 0x3000 && c <= 0x303f)) cjk++;
    }
    const nonCjk = text.length - cjk;
    return Math.floor(nonCjk / 4 + cjk / 1.8) + 1;
  };

  const inputTokens = useMemo(() => estimateInputTokens(input), [input]);

  const formatTokens = (n: number) => {
    if (n >= 1000) return (n / 1000).toFixed(1) + "K";
    return String(n);
  };

  const isAuthError = (msg: string) =>
    /api.?key|unauthorized|authentication|401|403|invalid.*key/i.test(msg);

  const loadHistory = () => {
    if (!currentSessionId) return;
    if (skipHistory.current) {
      skipHistory.current = false;
      return;
    }
    fetch(`/api/v1/chat/${currentSessionId}/history`)
      .then((r) => r.json())
      .then((d) => {
        setMessages(buildMessagesFromApi(d.messages || []));
      })
      .catch(() => {});
  };

  useEffect(() => {
    loadHistory();
  }, [currentSessionId]);

  useEffect(() => {
    messagesEnd.current?.scrollIntoView({ behavior: "smooth" });
  }, [messages, streamingSegments]);

  const messageQueueRef = useRef<string[]>([]);

  const doStream = async (text: string, sid: string) => {
    setStreamingSegments([]);
    setStreaming(true);

    const controller = new AbortController();
    abortRef.current = controller;

    const showError = (msg: string) => {
      setMessages((prev) => [
        ...prev,
        { id: `err_${Date.now()}`, role: "system", segments: [{ type: "text", content: msg }], isError: true },
      ]);
    };

    try {
      const resp = await fetch(`/api/v1/chat/${sid}/stream`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ content: text }),
        signal: controller.signal,
      });

      if (!resp.ok) {
        showError(`HTTP ${resp.status}: ${resp.statusText}`);
        return;
      }

      const reader = resp.body?.getReader();
      if (!reader) return;

      const decoder = new TextDecoder();

      while (true) {
        const { done, value } = await reader.read();
        if (done) break;

        const chunk = decoder.decode(value, { stream: true });
        const lines = chunk.split("\n");

        for (const line of lines) {
          if (!line.startsWith("data: ")) continue;
          const data = line.slice(6);
          if (data === "[DONE]") {
            setStreamingSegments((prev) => {
              if (prev.length > 0) {
                const usage = pendingUsageRef.current;
                pendingUsageRef.current = undefined;
                setMessages((msgs) => [
                  ...msgs,
                  {
                    id: `assistant_${Date.now()}`,
                    role: "assistant",
                    segments: prev,
                    tokenUsage: usage,
                  },
                ]);
              }
              return [];
            });
            break;
          }
          try {
            const parsed = JSON.parse(data);

            switch (parsed.type) {
              case "content_chunk":
                if (parsed.text) {
                  setStreamingSegments((prev) => {
                    const last = prev[prev.length - 1];
                    if (last && last.type === "text") {
                      const updated = [...prev];
                      updated[updated.length - 1] = {
                        ...last,
                        content: (last as TextSeg).content + parsed.text,
                      };
                      return updated;
                    }
                    return [...prev, { type: "text", content: parsed.text } as TextSeg];
                  });
                }
                break;

              case "tool_call_start":
                setStreamingSegments((prev) => [
                  ...prev,
                  {
                    type: "tool_call",
                    tool_name: parsed.tool_name,
                    args: parsed.args,
                    result: "",
                  } as ToolCallSeg,
                ]);
                break;

              case "tool_call_end":
                setStreamingSegments((prev) =>
                  prev.map((seg) =>
                    seg.type === "tool_call" &&
                    (seg as ToolCallSeg).tool_name === parsed.tool_name &&
                    !(seg as ToolCallSeg).result
                      ? { ...seg, result: parsed.result } as ToolCallSeg
                      : seg
                  )
                );
                break;

              case "done":
                if (parsed.prompt_tokens != null || parsed.usage) {
                  pendingUsageRef.current = {
                    prompt_tokens: parsed.prompt_tokens ?? parsed.usage?.prompt_tokens ?? 0,
                    completion_tokens: parsed.usage?.completion_tokens ?? 0,
                    total_tokens: parsed.usage?.total_tokens ?? 0,
                  };
                }
                setStreamingSegments((prev) => {
                  if (prev.length > 0) {
                    const usage = pendingUsageRef.current;
                    pendingUsageRef.current = undefined;
                    setMessages((msgs) => [
                      ...msgs,
                      {
                        id: `assistant_${Date.now()}`,
                        role: "assistant",
                        segments: prev,
                        tokenUsage: usage,
                      },
                    ]);
                  }
                  return [];
                });
                break;

              case "error":
                const errMsg = isAuthError(parsed.message)
                  ? t.error.auth
                  : parsed.message;
                showError(errMsg);
                break;
            }
          } catch {}
        }
      }
    } catch (err: any) {
      if (err.name === "AbortError") {
        setStreamingSegments((prev) => {
          if (prev.length > 0) {
            setMessages((msgs) => [
              ...msgs,
              {
                id: `assistant_${Date.now()}`,
                role: "assistant",
                segments: [...prev, { type: "text", content: `\n\n*(${t.chat.truncated})*` } as TextSeg],
              },
            ]);
          }
          return [];
        });
        return;
      }
      if (err.message?.includes("Failed to fetch") || err.name === "TypeError") {
        showError(t.error.network);
      } else {
        showError(`${t.error.unknown}: ${err.message}`);
      }
    } finally {
      setStreaming(false);
      abortRef.current = null;
      flushQueue();
    }
  };

  const flushQueue = () => {
    if (messageQueueRef.current.length > 0) {
      const next = messageQueueRef.current.shift()!;
      doStream(next, currentSessionId!);
    }
  };

  const send = async () => {
    const text = input.trim();
    if (!text) return;

    let sid = currentSessionId;
    if (!sid) {
      try {
        const res = await fetch("/api/v1/sessions", { method: "POST" });
        const data = await res.json();
        sid = data.id;
        skipHistory.current = true;
        setCurrentSessionId(sid);
        const updated = await fetch("/api/v1/sessions").then((r) => r.json());
        setSessions(updated.sessions || []);
      } catch {
        return;
      }
    }

    setInput("");

    const userMsg: Message = {
      id: `local_${Date.now()}`,
      role: "user",
      segments: [{ type: "text", content: text }],
    };
    setMessages((prev) => [...prev, userMsg]);

    if (streaming) {
      messageQueueRef.current.push(text);
      setMessages((prev) => [
        ...prev,
        {
          id: `hint_${Date.now()}`,
          role: "system",
          segments: [{ type: "text", content: t.chat.queued }],
        },
      ]);
      return;
    }

    doStream(text, sid);
  };

  const stop = () => {
    messageQueueRef.current = [];
    abortRef.current?.abort();
  };

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.nativeEvent.isComposing) return;
    if (e.key === "Enter" && !e.shiftKey) {
      e.preventDefault();
      send();
    }
  };

  const isLastStreamingText = (idx: number, segs: Segment[]) => {
    for (let i = idx + 1; i < segs.length; i++) {
      if (segs[i].type === "text") return false;
    }
    return true;
  };

  return (
    <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
      <div style={{ flex: 1, overflow: "auto", padding: "16px 24px" }}>
        {messages.length === 0 && !streaming && (
          <div
            style={{
              display: "flex", alignItems: "center", justifyContent: "center",
              height: "100%", color: "var(--text-muted)", fontSize: "14px",
            }}
          >
            {t.chat.empty}
          </div>
        )}

        {messages.map((msg) => (
          <MessageBlock key={msg.id} msg={msg} />
        ))}

        {streaming && streamingSegments.length === 0 && (
          <ThinkingBubble text={t.chat.thinking} />
        )}

        {streaming && streamingSegments.length > 0 && (
          <StreamingBlock segments={streamingSegments} isLastText={isLastStreamingText} />
        )}

        <div ref={messagesEnd} />
      </div>

        <InputBar
          input={input}
          setInput={setInput}
          streaming={streaming}
          onSend={send}
          onStop={stop}
          placeholder={t.chat.placeholder}
          sendLabel={t.chat.send}
          stopLabel={t.chat.stop}
          onKeyDown={handleKeyDown}
          inputTokens={inputTokens}
          formatTokens={formatTokens}
        />
    </div>
  );
}


function MessageBlock({ msg }: { msg: Message }) {
  if (msg.segments.length === 0) return null;
  const isUser = msg.role === "user";
  const isSystem = msg.role === "system";

  return (
    <div
      style={{
        marginBottom: "16px",
        display: "flex",
        justifyContent: isUser ? "flex-end" : isSystem ? "center" : "flex-start",
        flexDirection: "column",
        alignItems: isUser ? "flex-end" : isSystem ? "center" : "flex-start",
      }}
    >
      {msg.segments.map((seg, i) => (
        <SegmentView key={i} seg={seg} role={msg.role} isError={msg.isError} />
      ))}
      {!isUser && !isSystem && msg.tokenUsage && msg.tokenUsage.total_tokens > 0 && (
        <div style={{
          fontSize: "11px", color: "var(--text-muted)",
          marginTop: "2px", paddingLeft: "4px",
        }}>
          {msg.tokenUsage.prompt_tokens > 0 && (
            <span>{msg.tokenUsage.prompt_tokens.toLocaleString()} prompt</span>
          )}
          {msg.tokenUsage.completion_tokens > 0 && (
            <span>
              {msg.tokenUsage.prompt_tokens > 0 ? " + " : ""}
              {msg.tokenUsage.completion_tokens.toLocaleString()} completion
            </span>
          )}
          {msg.tokenUsage.total_tokens > 0 && (
            <span> = {msg.tokenUsage.total_tokens.toLocaleString()} total</span>
          )}
        </div>
      )}
    </div>
  );
}


function StreamingBlock({
  segments,
  isLastText,
}: {
  segments: Segment[];
  isLastText: (idx: number, segs: Segment[]) => boolean;
}) {
  return (
    <div
      style={{
        marginBottom: "16px",
        display: "flex",
        flexDirection: "column",
        alignItems: "flex-start",
      }}
    >
      {segments.map((seg, i) => {
        const isStreamingText = seg.type === "text" && isLastText(i, segments);
        return (
          <div key={i}>
            {seg.type === "text" ? (
              <TextBubble
                content={(seg as TextSeg).content}
                showCursor={isStreamingText}
              />
            ) : (
              <ToolCallBubble tc={seg as ToolCallSeg} />
            )}
          </div>
        );
      })}
    </div>
  );
}

function ThinkingBubble({ text }: { text: string }) {
  return (
    <div
      style={{
        marginBottom: "16px",
        display: "flex",
        justifyContent: "flex-start",
      }}
    >
      <div
        style={{
          padding: "10px 16px",
          borderRadius: "12px",
          fontSize: "14px",
          background: "var(--assistant-msg-bg)",
          color: "var(--text-secondary)",
          borderTopLeftRadius: "4px",
          display: "flex",
          alignItems: "center",
          gap: "4px",
        }}
      >
        {text}
        <span style={{ display: "inline-flex", gap: "3px" }}>
          {[0, 1, 2].map((i) => (
            <span
              key={i}
              style={{
                width: "5px",
                height: "5px",
                borderRadius: "50%",
                background: "var(--accent)",
                animation: `thinking-dot 1.2s ${i * 0.2}s infinite ease-in-out`,
              }}
            />
          ))}
        </span>
      </div>
    </div>
  );
}


function SegmentView({
  seg,
  role,
  isError,
}: {
  seg: Segment;
  role: string;
  isError?: boolean;
}) {
  if (seg.type === "text") {
    return (
      <TextBubble
        content={(seg as TextSeg).content}
        isUser={role === "user"}
        isError={isError}
        isSystem={role === "system"}
      />
    );
  }
  return <ToolCallBubble tc={seg as ToolCallSeg} />;
}


function TextBubble({
  content,
  showCursor,
  isUser,
  isError,
  isSystem,
}: {
  content: string;
  showCursor?: boolean;
  isUser?: boolean;
  isError?: boolean;
  isSystem?: boolean;
}) {
  const isMarkdown = !isUser && !isError && !isSystem;
  return (
    <div
      className={isMarkdown ? "markdown" : ""}
      style={{
        maxWidth: "70%",
        padding: "10px 16px",
        marginBottom: "4px",
        borderRadius: "12px",
        fontSize: "14px",
        lineHeight: "1.6",
        whiteSpace: isMarkdown ? "normal" : "pre-wrap",
        wordBreak: "break-word",
        background: isError
          ? "rgba(239,68,68,0.1)"
          : isUser
          ? "var(--user-msg-bg)"
          : "var(--assistant-msg-bg)",
        color: isError
          ? "var(--error)"
          : isUser
          ? "var(--user-msg-text)"
          : "var(--text-primary)",
        border: isError ? "1px solid var(--error)" : "none",
        borderTopRightRadius: isUser ? "4px" : "12px",
        borderTopLeftRadius: isUser ? "12px" : "4px",
      }}
    >
      {isUser || isError || isSystem ? (
        content
      ) : (
        <ReactMarkdown
          remarkPlugins={[remarkGfm]}
          rehypePlugins={[rehypeHighlight]}
        >
          {content}
        </ReactMarkdown>
      )}
      {showCursor && (
        <span
          style={{
            display: "inline-block",
            width: "8px", height: "16px",
            background: "var(--accent)",
            marginLeft: "2px",
            animation: "blink 1s infinite",
            verticalAlign: "middle",
          }}
        />
      )}
    </div>
  );
}


function InputBar({
  input, setInput, streaming, onSend, onStop,
  placeholder, sendLabel, stopLabel, onKeyDown,
  inputTokens, formatTokens,
}: any) {
  return (
    <div
      style={{
        padding: "12px 24px",
        borderTop: "1px solid var(--border-color)",
        background: "var(--bg-primary)",
      }}
    >
      <div style={{ display: "flex", gap: "8px", alignItems: "flex-end" }}>
        <textarea
          value={input}
          onChange={(e: any) => setInput(e.target.value)}
          onKeyDown={onKeyDown}
          placeholder={placeholder}
          rows={2}
          style={{
            flex: 1, resize: "none", padding: "10px 14px",
            borderRadius: "8px", border: "1px solid var(--border-color)",
            background: "var(--bg-secondary)", fontSize: "14px", lineHeight: "1.5",
          }}
        />
        {streaming ? (
          <button
            onClick={onStop}
            style={{
              padding: "10px 18px", borderRadius: "8px",
              background: "var(--error)", color: "white",
              fontSize: "14px", fontWeight: 500,
            }}
          >
            {stopLabel}
          </button>
        ) : (
          <button
            onClick={onSend}
            disabled={!input.trim()}
            style={{
              padding: "10px 18px", borderRadius: "8px",
              background: input.trim() ? "var(--accent)" : "var(--bg-tertiary)",
              color: input.trim() ? "white" : "var(--text-muted)",
              fontSize: "14px", fontWeight: 500,
            }}
          >
            {sendLabel}
          </button>
        )}
      </div>
      {input.length > 0 && (
        <div style={{
          marginTop: "4px", fontSize: "11px",
          color: "var(--text-muted)", textAlign: "right",
        }}>
          ~{formatTokens(inputTokens)} tokens
        </div>
      )}
    </div>
  );
}


function ToolCallBubble({ tc }: { tc: { tool_name: string; args: string; result: string } }) {
  const [expanded, setExpanded] = useState(false);
  const isComplete = !!tc.result;
  const isError = tc.result?.startsWith("Error:");

  let argsPreview = "";
  try {
    const parsed = JSON.parse(tc.args);
    argsPreview = Object.entries(parsed)
      .map(([k, v]) => `${k}=${JSON.stringify(v)}`)
      .join(", ")
      .slice(0, 80);
    if (argsPreview.length >= 80) argsPreview += "...";
  } catch {
    argsPreview = tc.args.slice(0, 80);
  }

  return (
    <div
      onClick={() => setExpanded(!expanded)}
      style={{
        marginBottom: "4px",
        padding: "6px 10px",
        borderRadius: "6px",
        fontSize: "12px",
        cursor: "pointer",
        background: isError
          ? "rgba(239,68,68,0.08)"
          : "rgba(99,102,241,0.08)",
        border: `1px solid ${isError ? "rgba(239,68,68,0.3)" : "rgba(99,102,241,0.2)"}`,
        color: "var(--text-secondary)",
        minWidth: "200px",
      }}
    >
      <div style={{ display: "flex", alignItems: "center", gap: "6px" }}>
        <span>{isComplete ? "🔧" : "⏳"}</span>
        <span style={{ fontWeight: 500 }}>{tc.tool_name}</span>
        <span style={{ opacity: 0.7 }}>
          {isComplete ? argsPreview : "executing..."}
        </span>
        <span style={{ marginLeft: "auto" }}>{expanded ? "▾" : "▸"}</span>
      </div>
      {expanded && (
        <div style={{ marginTop: "6px" }}>
          <div style={{ marginBottom: "4px", color: "var(--text-muted)" }}>
            args: {argsPreview}
          </div>
          {isComplete && (
            <div
              style={{
                padding: "6px 8px", borderRadius: "4px",
                background: "var(--bg-tertiary)", fontSize: "11px",
                whiteSpace: "pre-wrap", wordBreak: "break-word",
                maxHeight: "120px", overflow: "auto",
                color: isError ? "var(--error)" : "var(--text-secondary)",
              }}
            >
              {tc.result}
            </div>
          )}
        </div>
      )}
    </div>
  );
}

import { useEffect, useRef, useState } from "react";
import { useLocale } from "../../hooks/useLocale";
import { useStore } from "../../store";

export function ModelSelector() {
  const { t } = useLocale();
  const { models, currentModelId, setModels, setCurrentModelId } = useStore();
  const [open, setOpen] = useState(false);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState("");
  const ref = useRef<HTMLDivElement>(null);

  useEffect(() => {
    setLoading(true);
    setError("");
    Promise.all([
      fetch("/api/v1/llm/models").then((r) => r.json()),
      fetch("/api/v1/llm/current").then((r) => r.json()),
    ])
      .then(([modelsData, currentData]) => {
        setModels(modelsData.models || []);
        setCurrentModelId(currentData.id || "");
      })
      .catch(() => setError(t.model.error))
      .finally(() => setLoading(false));
  }, [t]);

  useEffect(() => {
    function handleClick(e: MouseEvent) {
      if (ref.current && !ref.current.contains(e.target as Node)) {
        setOpen(false);
      }
    }
    document.addEventListener("mousedown", handleClick);
    return () => document.removeEventListener("mousedown", handleClick);
  }, []);

  const switchModel = async (id: string) => {
    await fetch("/api/v1/llm/switch", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ model_id: id }),
    });
    setCurrentModelId(id);
    setOpen(false);
  };

  const current = models.find((m: any) => m.id === currentModelId);

  return (
    <div ref={ref} style={{ position: "relative" }}>
      <button
        onClick={() => setOpen(!open)}
        style={{
          padding: "6px 12px",
          borderRadius: "6px",
          fontSize: "13px",
          border: "1px solid var(--border-color)",
          background: "var(--bg-secondary)",
          display: "flex",
          alignItems: "center",
          gap: "6px",
        }}
      >
        {t.model.selector}:{" "}
        {loading
          ? t.model.loading
          : error
          ? t.model.error
          : current?.name || currentModelId || "?"}
        <span style={{ fontSize: "10px" }}>▼</span>
      </button>
      {open && (
        <div
          style={{
            position: "absolute",
            top: "100%",
            right: 0,
            marginTop: "4px",
            background: "var(--bg-primary)",
          border: `1px solid ${error ? "var(--error)" : "var(--border-color)"}`,
            borderRadius: "8px",
            minWidth: "200px",
            boxShadow: "0 4px 12px rgba(0,0,0,0.15)",
            zIndex: 100,
          }}
        >
          {models.map((m: any) => (
            <button
              key={m.id}
              onClick={() => switchModel(m.id)}
              style={{
                display: "block",
                width: "100%",
                padding: "10px 14px",
                textAlign: "left",
                fontSize: "13px",
                background:
                  m.id === currentModelId ? "var(--bg-hover)" : "none",
                borderBottom: "1px solid var(--border-color)",
              }}
            >
              <div>{m.name}</div>
              <div style={{ fontSize: "11px", color: "var(--text-muted)" }}>
                {m.provider} · context {m.limit?.context?.toLocaleString() || "?"}
              </div>
            </button>
          ))}
        </div>
      )}
    </div>
  );
}

import { useLocale } from "../../hooks/useLocale";
import { useStore } from "../../store";

export function SessionList() {
  const { t } = useLocale();
  const {
    sessions,
    currentSessionId,
    setCurrentSessionId,
    setSessions,
  } = useStore();

  const createSession = async () => {
    try {
      const res = await fetch("/api/v1/sessions", { method: "POST" });
      const data = await res.json();
      const updated = await fetch("/api/v1/sessions").then((r) => r.json());
      setSessions(updated.sessions || []);
      setCurrentSessionId(data.id);
    } catch (e) {
      console.error(e);
    }
  };

  const deleteSession = async (id: string, e: React.MouseEvent) => {
    e.stopPropagation();
    await fetch(`/api/v1/sessions/${id}`, { method: "DELETE" });
    const updated = await fetch("/api/v1/sessions").then((r) => r.json());
    setSessions(updated.sessions || []);
    if (currentSessionId === id) {
      setCurrentSessionId(updated.sessions?.[0]?.id || "");
    }
  };

  return (
    <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
      <div
        style={{
          padding: "12px 16px",
          borderBottom: "1px solid var(--border-color)",
        }}
      >
        <button
          onClick={createSession}
          style={{
            width: "100%",
            padding: "8px 0",
            borderRadius: "6px",
            background: "var(--accent)",
            color: "white",
            fontSize: "14px",
            fontWeight: 500,
          }}
        >
          + {t.sidebar.newChat}
        </button>
      </div>
      <div style={{ flex: 1, overflow: "auto", padding: "8px" }}>
        {sessions.length === 0 && (
          <div
            style={{
              padding: "16px",
              color: "var(--text-muted)",
              fontSize: "13px",
              textAlign: "center",
            }}
          >
            {t.sidebar.noSessions}
          </div>
        )}
        {sessions.map((s: any) => (
          <div
            key={s.id}
            onClick={() => setCurrentSessionId(s.id)}
            style={{
              padding: "10px 12px",
              marginBottom: "2px",
              borderRadius: "6px",
              cursor: "pointer",
              fontSize: "13px",
              background:
                s.id === currentSessionId
                  ? "var(--bg-hover)"
                  : "transparent",
              display: "flex",
              justifyContent: "space-between",
              alignItems: "center",
            }}
          >
            <span
              style={{
                overflow: "hidden",
                textOverflow: "ellipsis",
                whiteSpace: "nowrap",
                flex: 1,
              }}
            >
              {s.title}
            </span>
            <button
              onClick={(e) => deleteSession(s.id, e)}
              style={{
                marginLeft: "4px",
                padding: "2px 6px",
                fontSize: "11px",
                borderRadius: "4px",
                color: "var(--text-muted)",
                opacity: 0,
              }}
              onMouseEnter={(e) => (e.currentTarget.style.opacity = "1")}
              onMouseLeave={(e) => (e.currentTarget.style.opacity = "0")}
            >
              ✕
            </button>
          </div>
        ))}
      </div>
    </div>
  );
}

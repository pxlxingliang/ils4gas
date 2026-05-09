import { useLocale } from "../../hooks/useLocale";

export function LocaleToggle() {
  const { locale, toggleLocale } = useLocale();

  return (
    <button
      onClick={toggleLocale}
      style={{
        padding: "4px 8px",
        fontSize: "12px",
        borderRadius: "4px",
        border: "1px solid var(--border-color)",
        background: "var(--bg-secondary)",
        color: "var(--text-secondary)",
      }}
    >
      {locale === "zh-CN" ? "EN" : "中"}
    </button>
  );
}

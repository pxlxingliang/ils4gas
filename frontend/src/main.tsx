import React from "react";
import ReactDOM from "react-dom/client";
import { App } from "./App";
import { LocaleProvider } from "./hooks/useLocale";
import { initTheme } from "./store";
import "./styles/theme.css";
import "./styles/global.css";

initTheme();

ReactDOM.createRoot(document.getElementById("root")!).render(
  <React.StrictMode>
    <LocaleProvider>
      <App />
    </LocaleProvider>
  </React.StrictMode>
);

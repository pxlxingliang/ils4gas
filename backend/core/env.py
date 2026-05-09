import json
import os
from pathlib import Path
from typing import Dict


_DEFAULTS: Dict[str, str] = {
    "ILS4GAS_WEB_HOST": "0.0.0.0",
    "ILS4GAS_WEB_PORT": "8789",
    "ILS4GAS_MCP_TRANSPORT": "stdio",
    "ILS4GAS_MCP_HOST": "localhost",
    "ILS4GAS_MCP_PORT": "50001",
    "ILS4GAS_TRAIN_EB_PATH": "",
    "ILS4GAS_QM_FEATURE_PATH": "",
    "ILS4GAS_MOLECULE_GEN_SCRIPT": "",
    "ILS4GAS_EB_PREDICT_SCRIPT": "/personal/test/dwl/ils4gas-models/Model/Property_pred/Eb_predict.py",
}


class EnvManager:
    _initialized = False

    @classmethod
    def init(cls) -> None:
        if cls._initialized:
            return

        env_file = Path("~/.ils4gas/env.json").expanduser()
        existing: Dict[str, str] = {}
        if env_file.exists():
            try:
                existing = json.loads(env_file.read_text())
            except (json.JSONDecodeError, IOError):
                pass

        merged = dict(_DEFAULTS)
        merged.update(existing)

        changed = False
        for k, v in merged.items():
            if k not in existing:
                changed = True
            # Only set os.environ if not already explicitly set externally
            if k not in os.environ:
                os.environ[k] = v

        if changed or not env_file.exists():
            env_file.parent.mkdir(parents=True, exist_ok=True)
            env_file.write_text(json.dumps(merged, indent=2, ensure_ascii=False) + "\n")

        cls._initialized = True

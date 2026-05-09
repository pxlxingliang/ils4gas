# ILS4GAS

**AI Agent for Ionic Liquid Solubility of GAS** — a specialized framework for material science computations, with multiple user interfaces (Web UI / TUI / MCP server) for flexible access.

## Features

- **Crystal & Molecule Structure Generation** — Generate various crystal structures and molecules
- **Property Database Query** — Search ionic liquid and gas solubility properties
- **DFT Feature Retrieval** — Access DFT-calculated molecular features
- **Novel Ion Generation** — Generate novel ionic liquid molecules using machine learning
- **Binding Energy Prediction** — Predict binding energy between molecular structures

## Quick Start

### 1. Clone & Build Frontend

```bash
git clone https://github.com/pxlxingliang/ils4gas.git
cd ils4gas
cd frontend && npm install && npm run build && cd ..
```

### 2. Install

```bash
pip install ".[materials]"
```

Optional extras:

```bash
pip install ".[materials,tools]"   # includes openbabel + rdkit for SMILES conversion
```

### 3. Configure Environment (~/.ils4gas/env.json)

Environment variables are managed in `~/.ils4gas/env.json`. This file is created automatically on first run with default values.

Edit or create `~/.ils4gas/env.json`:

```json
{
  "ILS4GAS_TRAIN_EB_PATH": "/path/to/properties.pkl",
  "ILS4GAS_QM_FEATURE_PATH": "/path/to/dft_features.pkl",
  "ILS4GAS_MOLECULE_GEN_SCRIPT": "/path/to/generation_script.sh",
  "ILS4GAS_EB_PREDICT_SCRIPT": "/personal/test/dwl/ils4gas-models/Model/Property_pred/Eb_predict.py"
}
```

See [Environment Variables](#environment-variables) for the complete list.

### 4. Configure LLM (~/.ils4gas/config.json)

Copy the example config and set your API keys:

```bash
mkdir -p ~/.ils4gas
cp config.example.json ~/.ils4gas/config.json
# Edit ~/.ils4gas/config.json to add your API keys
```

Config example:

```json
{
  "currentModel": "volcengine/deepseek-v3.2",
  "providers": {
    "volcengine": {
      "name": "Volcano Engine",
      "type": "openai-compatible",
      "options": {
        "baseURL": "https://ark.cn-beijing.volces.com/api/coding/v3",
        "apiKey": "${ENV:VOLC_API_KEY}"
      },
      "models": {
        "deepseek-v3.2": {
          "name": "DeepSeek V3.2",
          "limit": { "context": 128000, "output": 4096 }
        }
      }
    }
  },
  "mcp_servers": [],
  "server": { "host": "0.0.0.0", "port": 8789 }
}
```

**Key points:**
- Environment variables in `${ENV:VAR_NAME}` format are resolved at runtime
- You can add multiple providers and switch between them at runtime
- `currentModel` sets the default model on startup

### 5. Choose Your Interface

#### TUI Mode (default)

```bash
ils4gas
```

#### Web UI Mode

```bash
ils4gas web
```

Open **http://localhost:8789** in your browser.

#### MCP Server Mode

Expose tools for other agents (Claude Code, etc.):

```bash
ils4gas mcp --transport stdio                        # stdio (for Claude Code)
ils4gas mcp --transport sse --port 50001             # SSE (network clients)
ils4gas mcp --transport streamable-http --port 50001 # HTTP streaming
```

## Environment Variables

All environment variables are managed via `~/.ils4gas/env.json`:

| Variable | Purpose | Default | Required For |
|----------|---------|---------|--------------|
| `ILS4GAS_WEB_HOST` | Web server host | `0.0.0.0` | Web UI |
| `ILS4GAS_WEB_PORT` | Web server port | `8789` | Web UI |
| `ILS4GAS_MCP_TRANSPORT` | MCP server transport | `stdio` | MCP server |
| `ILS4GAS_MCP_HOST` | MCP server host | `localhost` | MCP server |
| `ILS4GAS_MCP_PORT` | MCP server port | `50001` | MCP server |
| `ILS4GAS_TRAIN_EB_PATH` | Properties dataset path | — | Property search |
| `ILS4GAS_QM_FEATURE_PATH` | DFT features dataset path | — | DFT feature search |
| `ILS4GAS_MOLECULE_GEN_SCRIPT` | Ion generation script path | — | Novel ion generation |
| `ILS4GAS_EB_PREDICT_SCRIPT` | Binding energy prediction script | See note | Binding energy prediction |

**Note:** `ILS4GAS_EB_PREDICT_SCRIPT` defaults to `/personal/test/dwl/ils4gas-models/Model/Property_pred/Eb_predict.py`

## For Developers

See the `docs/` directory for developer documentation:

- **[docs/MCP_TOOLS.md](docs/MCP_TOOLS.md)** — How to add custom MCP tools
- **[docs/DESIGN.md](docs/DESIGN.md)** — Architecture and design decisions
- **[docs/STRUCTURE.md](docs/STRUCTURE.md)** — Codebase structure guide

## Dependencies

- Python >= 3.10
- fastapi + uvicorn (web server)
- openai (LLM client)
- mcp >= 1.9.0 (Model Context Protocol)
- ase (structure generation — optional, `[materials]` extra)
- Node.js >= 18 (frontend build only)

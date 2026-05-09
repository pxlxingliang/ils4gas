# Adding MCP Tools

MCP tools are defined as decorated Python functions in auto-discovered modules.

## Quick Start

1. Create a `.py` file under `backend/mcp_server/modules/` (must not start with `_`):

```
backend/mcp_server/modules/
├── materials.py    # existing
└── my_tools.py     # new — add your tools here
```

2. Write your tools:

```python
from backend.mcp_server import mcp


@mcp.tool()
def greet(name: str) -> str:
    """Greet someone by name.

    Args:
        name: The name to greet.

    Returns:
        A greeting string.
    """
    return f"Hello, {name}!"


@mcp.tool()
def add_numbers(a: float, b: float) -> float:
    """Add two numbers together.

    Args:
        a: First number.
        b: Second number.

    Returns:
        The sum of a and b.
    """
    return a + b
```

3. Restart the MCP server — tools are auto-discovered on startup.

## Rules

| Rule | Detail |
|------|--------|
| File location | `backend/mcp_server/modules/<name>.py` |
| File naming | Must not start with `_` |
| Decorator | Every tool function must use `@mcp.tool()` |
| Import | Always `from backend.mcp_server import mcp` |
| Type hints | Required on all parameters and return value |
| Docstring | Required — becomes the tool description visible to LLM clients |
| Parameters | Use `Literal["a", "b", "c"]` to restrict choices to enumerated values |

## Parameter Types

```python
from typing import Literal, Optional, List


@mcp.tool()
def example(
    required_str: str,                         # required string
    optional_int: int = 10,                    # optional with default
    choice: Literal["red", "green"] = "red",   # enumerated options
    flag: bool = False,                         # boolean flag
    items: List[str] = None,                   # list parameter
) -> dict:
    """Demonstrate all parameter types.

    Args:
        required_str: A required string parameter.
        optional_int: An optional integer (default 10).
        choice: Must be 'red' or 'green'.
        flag: A boolean flag.
        items: A list of strings.

    Returns:
        A dict summary.
    """
    return {"str": required_str, "int": optional_int, "choice": choice}
```

## Utility Functions

Place pure Python helpers (no `@mcp.tool()`) in `backend/mcp_server/util/`:

```
backend/mcp_server/util/
└── my_utils.py      # import from: backend.mcp_server.util.my_utils
```

Then import them in your tool module:

```python
from backend.mcp_server import mcp
from backend.mcp_server.util.my_utils import helper_function


@mcp.tool()
def my_tool(x: int) -> str:
    """Uses a helper."""
    return helper_function(x)
```

## Existing Utility Functions

`backend/mcp_server/util/comm.py` provides:

| Function | Purpose |
|----------|---------|
| `generate_work_path()` | Create a timestamped work directory |
| `run_command(cmd)` | Run a shell command with live output |
| `get_physical_cores()` | Count physical CPU cores |
| `xyz_to_smiles(path)` | Convert XYZ file to SMILES (needs openbabel) |
| `to_canonical_smiles(smiles)` | Normalize SMILES string (needs rdkit) |

## Testing Your Tools

Start the MCP server:

```bash
ils4gas mcp --transport stdio
```

Test with an MCP client (e.g., Claude Code) or via the MCP inspector:

```bash
npx @modelcontextprotocol/inspector python -m backend.mcp_server.server
```

## Environment Variables

Tools that need external data can read from environment variables:

```python
import os


@mcp.tool()
def lookup_data(key: str) -> str:
    """Look up data from a configured source.

    Args:
        key: The lookup key.

    Returns:
        The value for the key.
    """
    data_path = os.environ.get("ILS4GAS_DATA_PATH")
    if not data_path:
        raise ValueError("ILS4GAS_DATA_PATH not set")
    # ... load and query data
```

All ILS4GAS env vars use the `ILS4GAS_` prefix. Existing ones:

| Variable | Purpose |
|----------|---------|
| `ILS4GAS_TRANSPORT` | Transport protocol (stdio / sse / streamable-http) |
| `ILS4GAS_HOST` | Server host (default: localhost) |
| `ILS4GAS_PORT` | Server port (default: 50001) |
| `ILS4GAS_TRAIN_EB_PATH` | Path to properties dataset (.pkl) |
| `ILS4GAS_QM_FEATURE_PATH` | Path to DFT features dataset (.pkl) |
| `ILS4GAS_MOLECULE_GEN_SCRIPT` | Path to ion generation script |
| `ILS4GAS_EB_PREDICT_SCRIPT` | Path to binding energy prediction script |

## Built-in Material Science Tools

### `generate_bulk_structure`
Generate crystal structures via ASE (sc, fcc, bcc, hcp, diamond, zincblende, rocksalt)

### `generate_molecule_structure`
Generate molecule structures from ASE g2 collection or single atoms

### `search_properties_from_database`
Query ionic liquid-gas solubility properties (solvation free energy, binding energy, etc.)

### `search_dft_feature`
Retrieve DFT-calculated molecular features (ESP, polarizability, dipole moment, etc.)

### `generate_novel_ions`
Generate novel ionic liquid molecules using ML with FlashFormer

### `predict_binding_energy`
Predict binding energy between two molecular structures using pre-trained SchNet model

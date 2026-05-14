from typing import Dict, List, Optional

try:
    import tiktoken

    _HAS_TIKTOKEN = True
except ImportError:
    _HAS_TIKTOKEN = False

_MODEL_TO_ENCODING: Dict[str, str] = {
    "gpt-4o": "o200k_base",
    "gpt-4o-mini": "o200k_base",
    "gpt-4-turbo": "cl100k_base",
    "gpt-4": "cl100k_base",
    "gpt-3.5-turbo": "cl100k_base",
    "text-embedding-ada-002": "cl100k_base",
    "text-embedding-3-small": "cl100k_base",
    "text-embedding-3-large": "cl100k_base",
}

FALLBACK_ENCODING = "cl100k_base"


def _get_encoding(model_name: str):
    if not _HAS_TIKTOKEN:
        return None

    encoding_name = _MODEL_TO_ENCODING.get(model_name)
    if encoding_name is None:
        if model_name.startswith("gpt-4"):
            encoding_name = "o200k_base" if "o" in model_name else "cl100k_base"
        elif model_name.startswith("gpt-3.5"):
            encoding_name = "cl100k_base"
        else:
            encoding_name = FALLBACK_ENCODING

    try:
        return tiktoken.get_encoding(encoding_name)
    except Exception:
        return tiktoken.get_encoding(FALLBACK_ENCODING)


def _char_estimate(text: str) -> int:
    chars = len(text)
    cjk = sum(1 for c in text if "\u4e00" <= c <= "\u9fff" or "\u3000" <= c <= "\u303f")
    non_cjk = chars - cjk
    return int(non_cjk / 4.0 + cjk / 1.8) + 1


def _count_text_tokens(encoding, text: str) -> int:
    if encoding is not None:
        return len(encoding.encode(text))
    return _char_estimate(text)


def count_message_tokens(message: Dict, model_name: Optional[str] = None) -> int:
    encoding = _get_encoding(model_name) if model_name else None

    tokens = 3
    for key in ("role", "content", "name", "tool_call_id", "tool_calls"):
        if key in message:
            tokens += 1

    content = message.get("content")
    if isinstance(content, str):
        tokens += _count_text_tokens(encoding, content)
    elif isinstance(content, list):
        for part in content:
            if isinstance(part, dict):
                part_type = part.get("type", "")
                if part_type in ("text", "input_text"):
                    tokens += _count_text_tokens(encoding, part.get("text", ""))
                elif part_type == "image_url":
                    tokens += 85

    tool_calls = message.get("tool_calls")
    if isinstance(tool_calls, list):
        for tc in tool_calls:
            tokens += 4
            func = tc.get("function", {})
            tokens += _count_text_tokens(encoding, func.get("name", ""))
            tokens += _count_text_tokens(encoding, func.get("arguments", ""))

    return tokens


def count_tokens(
    messages: List[Dict], model_name: Optional[str] = None
) -> int:
    total = 3
    for msg in messages:
        total += count_message_tokens(msg, model_name)
    return total


def estimate_text_tokens(text: str) -> int:
    return _char_estimate(text)

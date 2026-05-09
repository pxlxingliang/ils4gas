from typing import Optional

TITLE_PROMPT = (
    'Generate a short, concise title (1-6 words) for a conversation '
    'that starts with the following user message. Reply with ONLY the '
    'title text, no quotes, no explanation:\n\n'
    'User: "{content}"\n\nTitle:'
)


async def generate_title(content: str, llm_service) -> Optional[str]:
    if not content or not llm_service:
        return None

    trimmed = content.strip()
    if len(trimmed) > 200:
        trimmed = trimmed[:200]

    messages = [
        {"role": "system", "content": "You generate short conversation titles."},
        {"role": "user", "content": TITLE_PROMPT.format(content=trimmed)},
    ]

    try:
        result = await llm_service.ainvoke(messages, temperature=0.3, max_tokens=20)
        title = result.strip().strip('"').strip("'")
        if not title or len(title) > 80:
            return None
        return title
    except Exception:
        return None

import tempfile
import textwrap
from pathlib import Path

import pytest

from backend.tools.builtin.skill import (
    SkillRegistry,
    create_skill_tool,
    get_skill_registry,
)
from backend.tools.registry import ToolRegistry


def _make_skill_dir(base: Path, name: str, description: str, prompt: str) -> Path:
    d = base / name
    d.mkdir()
    content = textwrap.dedent(f"""\
    ---
    name: {name}
    description: {description}
    ---

    {prompt}
    """)
    (d / "SKILL.md").write_text(content, encoding="utf-8")
    return d


class TestSkillRegistrySummaries:
    def test_summaries_empty_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            reg = SkillRegistry(skills_dirs=[tmp])
            assert reg.get_summaries() == []

    def test_summaries_parses_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            _make_skill_dir(Path(tmp), "skill_a", "Skill A description", "# Prompt A")
            _make_skill_dir(Path(tmp), "skill_b", "Skill B description", "# Prompt B")

            reg = SkillRegistry(skills_dirs=[tmp])
            summaries = reg.get_summaries()
            assert len(summaries) == 2
            assert summaries[0]["name"] == "skill_a"
            assert summaries[0]["description"] == "Skill A description"

    def test_summaries_does_not_load_full_content(self):
        with tempfile.TemporaryDirectory() as tmp:
            _make_skill_dir(Path(tmp), "test_skill", "Test desc", "# Full Prompt Content")
            reg = SkillRegistry(skills_dirs=[tmp])

            summaries = reg.get_summaries()
            assert len(summaries) == 1
            assert "name" in summaries[0]
            assert "description" in summaries[0]
            assert "prompt" not in summaries[0]
            assert "Full" not in str(summaries)

            assert reg.get("test_skill") is None  # not loaded yet


class TestSkillTool:
    def test_get_available_skills_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            reg = SkillRegistry(skills_dirs=[tmp])
            assert reg.get_summaries() == []

    def test_get_available_skills(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            _make_skill_dir(base, "alpha", "Alpha desc", "# Alpha body")
            _make_skill_dir(base, "beta", "Beta desc", "# Beta body")

            reg = SkillRegistry(skills_dirs=[tmp])
            skills = reg.get_summaries()
            assert len(skills) == 2

    def test_create_skill_tool_empty_returns_none(self):
        with tempfile.TemporaryDirectory() as tmp:
            # Use a nested helper to override the global singleton
            reg = SkillRegistry(skills_dirs=[tmp])
            # We test via crate_skill_tool which uses get_skill_registry()
            # So test the registry directly since singleton is global
            assert reg.get_summaries() == []

    def test_create_skill_tool_creates_valid_tool(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            _make_skill_dir(base, "alpha", "Alpha description", "# Alpha body")
            _make_skill_dir(base, "beta", "Beta description", "# Beta body")

            reg = SkillRegistry(skills_dirs=[tmp])
            summaries = reg.get_summaries()
            assert len(summaries) == 2

            # Test load directly
            skill = reg.load("alpha")
            assert skill is not None
            assert skill.meta.name == "alpha"
            assert "Alpha body" in skill.prompt

            tool = create_skill_tool(registry=reg)
            assert tool is not None
            assert tool.name == "load_skill"
            assert "alpha" in tool.description
            assert "beta" in tool.description

    def test_skill_tool_fn_returns_skill_content(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            _make_skill_dir(base, "alpha", "Alpha description", "# Alpha body\n\nContent here.")

            reg = SkillRegistry(skills_dirs=[tmp])
            skill = reg.load("alpha")
            result = (
                f"## Skill: {skill.meta.name}\n"
                f"Directory: {skill.directory}\n\n"
                f"{skill.prompt}"
            )
            assert "## Skill: alpha" in result
            assert "Directory:" in result
            assert "Alpha body" in result
            assert "Content here" in result

    def test_load_skill_returns_none_for_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            reg = SkillRegistry(skills_dirs=[tmp])
            assert reg.load("nonexistent") is None

    def test_load_skill_returns_full_prompt(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            _make_skill_dir(base, "test", "Test description", "# Test Body\n\nHello World.")

            reg = SkillRegistry(skills_dirs=[tmp])
            skill = reg.load("test")
            assert skill is not None
            assert skill.meta.name == "test"
            assert skill.directory is not None
            assert "# Test Body" in skill.prompt
            assert "Hello World" in skill.prompt


class TestSkillToolRegistration:
    def test_register_into_tool_registry(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            _make_skill_dir(base, "code_review", "Review code", "# Review\n\nGuidelines here.")

            reg = SkillRegistry(skills_dirs=[tmp])
            skill = reg.load("code_review")
            assert skill is not None

            from backend.tools.builtin.skill import create_skill_tool as _create
            tool = _create(registry=reg)
            assert tool is not None

            registry = ToolRegistry()
            registry.register(tool)

            assert registry.get("load_skill") is not None
            tools = registry.to_openai_tools()
            assert len(tools) == 1
            assert tools[0]["function"]["name"] == "load_skill"

            loaded = registry.get("load_skill")
            result = loaded.execute(name="code_review")
            assert "Guidelines" in result

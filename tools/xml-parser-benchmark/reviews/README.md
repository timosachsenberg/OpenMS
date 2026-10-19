# Reviewer evidence

`review-input.txt` is the exact 212,224-byte full input used by the three completed vendor reviews. Its SHA-256 is recorded in each run JSON. `*-findings.md` contains final answers; these are unverified reviewer assertions until checked against `REVIEW-DECISIONS.md` and measured results.

The saved run metadata records CLI versions, settings, input hashes, elapsed times and outcomes. Failed Vibe attempts have no attributed findings. The successful Vibe run used the installed default requested model alias with thinking off in a task-local profile; its actual serving model was not independently identified. Claude's reported model-usage key is retained. Kimi's serving model is not identified by the saved final output.

With the matching CLI versions and existing authentication, run `python3 run_reviewer.py claude` or `python3 run_reviewer.py kimi` from this directory. To recreate Vibe's successful profile layout:

```sh
mkdir -p .vibe/agents
cp vibe-reviewer.toml .vibe/agents/optimized-review.toml
python3 run_vibe_final.py
```

That runner uses a per-process agent search path and disables every tool; it does not require changing global configuration or trusting the project. The generic Vibe invocation in `run_reviewer.py` reproduces the earlier full high-thinking attempt, not the successful thinking-off run. `vibe-focused-input.txt` and its runner describe a narrower failed attempt, also excluded from findings.

Private reasoning output and Kimi stderr are deliberately excluded. Only completed assistant final text was extracted for the published findings. No model was permitted to run benchmarks or edit the implementations; proposed changes were implemented and tested separately.

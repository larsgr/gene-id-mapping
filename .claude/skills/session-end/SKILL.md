---
name: session-end
description: This skill should be used when the user says "session-end", "/session-end", "document this session", "end session", or when work appears complete and needs documentation. Use this skill to document AI session work and prepare for commit.
version: 1.0.0
---

# Session-End Documentation Workflow

This skill guides you through documenting an AI session in the gene-id-mapping project following the established three-tier workflow structure.

## When This Skill Applies

Use this skill when:
- User explicitly requests `/session-end` or "document this session"
- A work session is complete and needs documentation
- You've finished implementing features or fixes
- User says "end session" or similar

## Workflow Steps

Follow these steps in order:

### 1. Identify Roadmap Task

Read the README.md file to identify current Roadmap tasks marked with 🔄 (in progress).

Ask the user which Roadmap task this session relates to:
- List tasks from README.md with 🔄 status
- Let user select one or specify "Other"
- Determine the detailed log file: `docs/[task-name].md`
  - If task log doesn't exist yet, you'll create it
  - Use descriptive names (e.g., `comparison-script-refine-and-debug.md`)

### 2. Gather Session Information

Use AskUserQuestion to collect:

**Required information:**
- **Session title**: Short description (will be used in commit message and all documentation)
  - Example: "Add multi-file input support"
  - Example: "Fix CDS phase calculation bug"
  - Should be imperative form (not past tense)

- **What was accomplished**: Bullet points of actions taken
  - Commands executed
  - Files created or modified
  - Analysis performed
  - Tests run

**Optional information:**
- **Key findings or conclusions**: Important discoveries or results
- **Problems encountered**: Issues faced and how they were resolved
- **User's original prompt**: The initial request that started this session

### 3. Generate Documentation Updates

Create entries for both documentation files:

#### For docs/AI-usage.md (append to end):

```markdown
### (Claude Code) [Session Title] (YYYY-MM-DD)

**Prompt:**
\`\`\`
[User's original request or prompt]
\`\`\`

**What it did:**
* [Action 1: command or file operation]
* [Action 2: what was created/modified]
* [Action 3: analysis or testing performed]
* [Include any errors encountered and resolutions]

**Reflection:**
[Your assessment of what worked well and what could be improved]
[Note technical challenges, error recovery strategies, suggestions for future work]
```

#### For docs/[task-name].md (append new section):

```markdown
## [Session Title] (YYYY-MM-DD)

**What was done:**
* [Bullet point list of actions taken]
* [Commands executed]
* [Files created/modified]
* [Analysis performed]

**Findings:**
[Key conclusions, discoveries, or results from this work]

**Issues:**
[Problems encountered and how they were resolved]
[Leave this section out if no significant issues]
```

**Important notes:**
- Use today's date in YYYY-MM-DD format
- Session title should match across both files and the commit message
- The AI-usage.md entry should include your reflection on the work
- The detailed task log should be more focused on technical outcomes

### 4. Show Preview and Confirm

Display what will be added to both files:
1. Show the complete entry for docs/AI-usage.md
2. Show the complete section for docs/[task-name].md
3. Ask: "Does this documentation look correct? Should I proceed with updating the files?"

If user requests changes, revise and show preview again.

### 5. Update Documentation Files

Once confirmed:
- Append new entry to docs/AI-usage.md (at the end of file)
- Append new section to docs/[task-name].md (at the end of file)
- If task log file doesn't exist, create it with:
  ```markdown
  # [Task Name]

  This document tracks detailed progress on the [task description] task from the README.md Roadmap.

  ## [Session Title] (YYYY-MM-DD)
  [... rest of first session entry ...]
  ```

### 6. Suggest Commit Message and Show Checklist

Present commit preparation information:

```
✅ Documentation updated:
   - docs/AI-usage.md (new entry added)
   - docs/[task-name].md (new section added)

📋 Commit Checklist:
   [ ] All documentation complete
   [ ] Large files (>1MB) added to .gitignore
   [ ] All intended changes staged
   [ ] Commit message matches session name

Suggested commit message:
"[Session Title]"

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>
```

### 7. Offer to Create Commit

Ask: "Would you like me to stage the changes and create the commit now?"

If yes:
1. Run `git status` to show what will be committed
2. Stage relevant files (documentation and work files)
3. Create commit with suggested message including Co-Authored-By line
4. Run `git status` again to confirm

If no:
- Remind user of the suggested commit message
- Note that they can commit manually when ready

## Best Practices

- **Be thorough**: Capture enough detail that someone reading the log later can understand what was done and why
- **Be concise**: Don't include unnecessary verbosity or redundant information
- **Match naming**: Session title must be consistent across AI-usage.md, task log, and commit message
- **Reflect honestly**: In AI-usage.md reflection, note what worked well AND what could be improved
- **Date format**: Always use YYYY-MM-DD format for dates
- **Imperative mood**: Use imperative form for session titles (e.g., "Add feature" not "Added feature")

## Special Cases

### First Session for a New Task
If the detailed task log doesn't exist yet:
1. Ask user for a descriptive filename
2. Create the file with a header explaining the task
3. Add the first session section

### Multiple Tasks in One Session
If work touched multiple Roadmap tasks:
- Ask user which is the primary task
- Document in that task's log
- Optionally mention other affected tasks in the documentation

### Session Without Original Prompt
If user started with conversation or exploration:
- In the "Prompt" section, write a brief summary of what initiated the work
- Or use "[Exploratory session - no specific initial prompt]"

## Integration with Git Workflow

This skill integrates with the project's git commit workflow:
- Documentation updates should be committed together with code changes
- Commit message should match the session title
- Always include Co-Authored-By line for Claude
- Check git status before and after committing

## Error Handling

If something goes wrong:
- If files can't be read: Ask user to check file paths
- If log file doesn't exist: Offer to create it
- If user skips required information: Prompt again with explanation
- If preview is rejected: Revise based on feedback and show again

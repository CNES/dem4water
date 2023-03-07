# **dem4water** **Contributing guide**.

1. [Bug report](#bug-report)
2. [Contributing workflow](#contributing-workflow)
3. [Contribution license agreement](#contribution-license-agreement)
4. [Coding guide](#coding-guide)
5. [Pylint pre-commit validation](#pylint-pre-commit-validation)
6. [Merge request acceptation process](#merge-request-acceptation-process)

# Bug report

Any proven or suspected malfunction should be traced in a bug report, the latter being an issue in the dem4water github repository.

**Don't hesitate to do so: It is best to open a bug report and quickly resolve it than to let a problem remains in the project.**
**Notifying the potential bugs is the first way for contributing to a software.**

In the problem description, be as accurate as possible. Include:

- The procedure used to initialize the environment
- The incriminated command line or python function

# Contributing workflow

Any code modification requires a Merge Request. It is forbidden to push patches directly into master (this branch is protected).

It is recommended to open your Merge Request as soon as possible in order to inform the developers of your ongoing work.
Please add `WIP:` before your Merge Request title if your work is in progress: This prevents an accidental merge and informs the other developers of the unfinished state of your work.

The Merge Request shall have a short description of the proposed changes. If it is relative to an issue, you can signal it by adding `Closes xx` where xx is the reference number of the issue.

Likewise, if you work on a branch (which is recommended), prefix the branch's name by `xx-` in order to link it to the xx issue.

dem4water Classical workflow is :

- Check Licence and sign [Contributor Licence Agreement](#contribution-license-agreement) (Individual or Corporate)
- Create an issue (or begin from an existing one)
- Create a Merge Request from the issue: a MR is created accordingly with "WIP:", "Closes xx" and associated "xx-name-issue" branch
- Git add, commit and push from local working clone directory or from the forge directly
- Follow [Conventional commits](https://www.conventionalcommits.org/) specifications for commit messages
- Beware that quality pre-commit tools are installed in continuous integration with classical quality code tools.
- When finished, change your Merge Request name (erase "WIP:" in title ) and ask `@dem4water` to review the code (see below Merge request acceptation process)

# Contribution license agreement

dem4water requires that contributors sign out a [Contributor License
Agreement](https://en.wikipedia.org/wiki/Contributor_License_Agreement). The
purpose of this CLA is to ensure that the project has the necessary ownership or
grants of rights over all contributions to allow them to distribute under the
chosen license (Apache License Version 2.0)

To accept your contribution, we need you to complete, sign and email to _dem4water [at]
cnes [dot] fr_ an [Individual Contributor Licensing
Agreement](./docs/source/CLA/ICLA-dem4water.doc) (ICLA) form and a
[Corporate Contributor Licensing
Agreement](./docs/source/CLA/CCLA-dem4water.doc) (CCLA) form if you are
contributing on behalf of your company or another entity which retains copyright
for your contribution.

The copyright owner (or owner's agent) must be mentioned in headers of all modified source files and also added to the [AUTHORS.md
file](./AUTHORS.md).

# Merge request acceptation process

The Merge Request will be merged into master after being reviewed by dem4water steering committee (core committers) composed of:

- Santiago Pena Luque (CNES)

Only the members of this committee can merge into master.

The checklist of a Merge Request acceptance is the following:

- [ ] At least one code review has been done by members of the group above (who check among others the correct application of the rules listed in the [Coding guide](#coding-guide)).
- [ ] All comments of the reviewers has been dealt with and are closed
- [ ] The reviewers have signaled their approbation (thumb up)
- [ ] No reviewer is against the Merge Request (thumb down)

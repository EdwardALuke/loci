# Contributing

Contributions to Loci are welcome. Loci currently uses a [public Github
repository][1] to manage collaboration.

This document provides a basic outline of the contribution process and the
types of contributions needed at this stage. A more in-depth Developer Guide
is on the roadmap.

[1]: https://github.com/EdwardALuke/loci

## Reporting

Once you have identified some need, the first step is to communicate it to the
developers by opening a new issue on Github. We recommend using the "Labels"
menu on the right hand side of the `New Issue` page to tag it with the
appropriate category. The three primary categories of issue we expect to see
are:

1. **Bug**: Apparent errors in the code or things which don't function as
   expected.
2. **Documentation**: Suggested corrections, improvements, or additions to the
   documentation.
3. **Enhancement**: Requests for new features or functionality.

Ensure your issue is described in enough detail for others to follow. If you
believe you have discovered a bug please provide any error messages, the steps
necessary to recreate the error, and the behavior you expected to see. Please
ensure you respond to any follow-up questions from the developers on your
issue in a timely fashion. Try to keep each issue as narrow in scope as is
reasonable. If your issue touches on multiple aspects of the code base, we
prefer that it be split into multiple sub-issues.

This is the minimum requirement for a contribution, and you are welcome to
leave it in the hands of the developers from here. If you would like to address
it yourself, continue on to the next section.

## Pull Requests

You are welcome to work on any open issue, whether you or someone else reported
it. To avoid duplicated effort, it is a good idea to leave a comment on the
issue stating that you are volunteering to work on it.

The general process for contributing via Github is:

1. Fork the Loci repo and create your issue branch from `dev`. See [2].
2. Create a draft pull request which would merge your issue branch on your
   personal fork into the `dev` branch on the Loci repository.
3. Work on fixing the issue in your personal fork of Loci.
4. Periodically merge the original Loci `dev` branch into your issue branch to
   avoid diverging, resolving any merge conflicts as they come up.
5. Run tests (see [Testing](#testing)) to check whether anything was obviously
   broken.
6. When you believe your contribution is ready, mark the pull request "ready
   for review."
7. Address any feedback in code review promptly.

If accepted, pull requests will be "squashed" (collapsed into a single
monolithic commit) and merged into the `dev` branch. This makes changes easier
to revert if needed, and it means developers don't need to be concerned about
the number of commits they are making, so commit early and often. This also
underlines the importance of keeping your issues narrow in scope and touching
the minimum number of files necessary to resolve the issue.

[2]: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo

## Testing

Loci's testing infrastructure is a work in progress, but should be used as a
sanity check on potential contributions.

[TBA: How to run tests.]

## General Guidelines

- **Commit messages**: The title of your commit message should be <= 50
  characters and start with the number of the issue you're working on, which
  tells Github to link your commit to the Issue. We generally adopt the style
  `[#108]:`. This is then followed by a short description of what was done
  written in imperative present tense, e.g. `Rename X to Y` or
  `Add [feature] to [module]`. Optionally, you may add a detailed description
  following a blank line after the title.

- **Code style**: In lieu of a concrete style guide, try to keep your code formatted
  similarly to the surrounding code. Formal guidance will be included in the
  Developer Guide when it is finished.

- **Fixing-as-you-go**: We are aware that many parts of the code base would
  benefit from small changes, such as extra comments, typo fixes, and
  formatting consistency. Instead of fixing these as you come across them,
  please open up issues which document them as you go, and follow the process
  for fixing them with a new pull request. This helps maintain traceability and
  avoid merge conflicts.

## License

By contributing to Loci, you agree that your contributions will be licensed
under the LGPL-3.0 license as provided in [COPYING.LESSER](COPYING.LESSER).

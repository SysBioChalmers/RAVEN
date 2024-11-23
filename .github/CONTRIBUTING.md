## Contributor guidelines

Anybody is welcome to contribute to the development of the RAVEN Toolbox, but please abide by the following guidelines.

Each function should start with a commented section describing the function and explaining the parameters. Existing functions can clarify what style should be used.

### Bugfixes, New Features and Functions
* For any development, whether bug fixes, introducing new functions or new/updated features for existing functions: make a separate branch from `devel` and name the branch for instance after the function/feature you are fixing/developing. If you work on a fix, start the branch name with `fix/`, if you work on a feature, start the branch name with `feat/`. Examples: `fix/format_reactions` or `feat/new_algorithms`.
* Make commits to this branch while developing. Aim for backward compatibility, and try to avoid very new MATLAB functions when possible, to accommodate users with older MATLAB versions.
* When you are happy with your new function/feature, make a pull request to the `devel` branch. Also, see [Pull request](#pull-request) below.

### Semantic Commits
Use semantic commit messages to make it easier to show what you are aiming to do:
* `chore`: updating binaries, KEGG or MetaCyc database files, etc.
* `doc`: updating documentation (in `doc` folder) or explanatory comments in functions.
* `feat`: new feature added, e.g. new function introduced / new parameters / new algorithm / etc.
* `fix`: bugfix.
* `refactor`: see [code refactoring](https://en.wikipedia.org/wiki/Code_refactoring).
* `style`: minor format changes of functions (spaces, semi-colons, etc., no code change).

Examples:
```
feat: exportModel additional export to YAML
chore: update KEGG model to version 83.0
fix: optimizeProb parsing results from Gurobi
```
A more detailed explanation or comments can be left in the commit description.

### External Software Update
* Once the newer version for any external software (BLAST+, CD-HIT, DIAMOND, HMMER, MAFFT, WoLFPSORT) is available, identify the newest version which is simultaneously available for macOS, Unix/Linux and Windows
* Create a separate branch from the `devel` branch and name it e.g. `chore/updateBinaries`
* Commit the changes for each program and operating system separately, e.g. `chore: update HMMER (Win) to 3.2.1`
* As soon as all binaries for particular program are updated through commits, update the corresponding license file, if available and place it in e.g. `software/blast+` directory.
* Update the version number in `RAVENdir/software/version.txt` file.
* Do some testing to ensure that the new binaries are working correctly. Upon successful tests, create a Pull Request to the `devel` branch.

### Pull Requests
* No changes should be directly committed to the `main` or `devel` branches. Commits are made to side branches, after which pull requests are made for merging with `main` or `devel`.
* The person making the pull request and the one accepting the merge _cannot_ be the same person.
* Typically, wait a few days before merging, to allow for other developers to inspect the pull request.
* A merge with the mainbranch invokes a new release (see [versioning](#versioning)).

### Versioning
RAVEN Toolbox follows [semantic versioning](https://semver.org/). After each successful Pull Request (PR) to `main`, the `version.txt` file should be updated in a separate PR. When this PR is merged, a new [release](https://github.com/SysBioChalmers/RAVEN/releases) should be specified, with a list detailing the most significant changes since the previous release.

Ensure that the documentation is also updated (in `devel`) before creating a PR to `main`. It can be updated by running `updateDocumentation()`.

### Categorize Releases
* PRs pushed into `main` should be patch releases by default.
* Developers who push PRs to `devel` may propose target release categories (e.g. minor or patch) that can be commented on and discussed by reviewers and other developers.
* As a result of the above step, minor releases hopefully can be determined.
* To help the process of categorizing releases, some conditions are suggested qualification for minor releases:
    * with database updates (e.g. KEGG, MetaCyc)
    * script reusability (i.e. with parameter changes in RAVEN functions)
    * with binary file changes

## Make a New Release
* Start a PR from develop to main with summary of all PRs since the last release, titled "RAVEN 2.9.3" (replace 2.9.3 with the version of the next release)
* Once ready, locally merge develop into main by accepting all changes using following commands, and push to origin:
    * `git checkout main`
    * `git pull`
    * `git merge -Xtheirs develop -m "RAVEN 2.9.3"` (replace 2.9.3 with the version of the next release)
    * `printf "2.9.3" > version.txt` (replace 2.9.3 with the version of the next release)
    * `git stage version.txt`
    * `git commite --amend --no-edit`
    * `git push`
* Initiate a new release, with name and tag "v2.9.3" (replace 2.9.3 with the version of the next release)
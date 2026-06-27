# nrf1-proteomics

[![ci](https://github.com/br3ndonland/nrf1-proteomics/workflows/ci/badge.svg)](https://github.com/br3ndonland/nrf1-proteomics/actions/workflows/ci.yml)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

## Description

**Welcome!** This is a template repository for Python projects, engineered for use as a [GitHub template repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-repository-from-a-template). To use the template, click on "[Use this template](https://github.com/new?template_name=nrf1-proteomics&template_owner=br3ndonland)" or use the [GitHub CLI](https://cli.github.com/manual/) (`gh repo create your-repo-name --public --template br3ndonland/nrf1-proteomics`).

The template repo name can be replaced as shown in the [quickstart](#quickstart) section.

## Quickstart

[Install Hatch](https://hatch.pypa.io/latest/install/), rename the project, then install the project:

```sh
# set your names
repo_name="your-repo-name"
your_name="Your Name"
your_user="YourGitHubUsername"

# update repo for new names
git mv "src/nrf1_proteomics" "src/${repo_name//-/_}"
git grep -l "br3ndonland" | xargs sed -i "s|br3ndonland|$your_user|g"
git grep -l "Brendon Smith" | xargs sed -i "s|Brendon Smith|$your_name|g"
git grep -l "nrf1_proteomics" | xargs sed -i "s|nrf1_proteomics|${repo_name//-/_}|g"
git grep -l "nrf1-proteomics" | xargs sed -i "s|nrf1-proteomics|$repo_name|g"

# run tests to verify
hatch run coverage run
```

## Documentation

Documentation is built with [Zensical](https://zensical.org/).

The [Vercel project configuration](https://vercel.com/docs/project-configuration) in `vercel.json` can be used to deploy the docs to [Vercel](https://vercel.com/docs).

## Related

Another common approach, especially for Python, is to use [cookiecutter](https://github.com/cookiecutter/cookiecutter). In a cookiecutter repo, the developer adds template variables throughout, like `{{cookiecutter.repo_name}}`. When a user runs `cookiecutter` using the template repository, the template variables are replaced with the information the user provides. This repo is simple enough that I haven't needed to add cookiecutter yet.

[Copier](https://copier.readthedocs.io/en/stable/) and [PyScaffold](https://pyscaffold.org/en/stable/) are similar to cookiecutter, with some additional benefits.

# Publishing this repository

## 1. Create the repository (you, not me)

I have not created anything on GitHub and I do not need a token or any other credential — please don't
paste a personal access token into a chat. Everything here is a local git repository with one commit, so
you can publish it yourself in under a minute:

```bash
cd CAD-scV2E2G
gh repo create jamrute/CAD-scV2E2G --private --source=. --push      # GitHub CLI
# or, with the web UI: create an empty repo named CAD-scV2E2G, then
git remote add origin git@github.com:jamrute/CAD-scV2E2G.git
git push -u origin main
```

Keep it private until the paper is out, then flip it to public. The manuscript already cites
`https://github.com/jamrute/CAD-scV2E2G`, so the name should stay as is (and the stray space in the Code
Availability sentence needs removing).

If you would rather I push, the safe route is a fine-grained token scoped to this one repository with
contents write only, expiring within a day — but doing it yourself with `gh` is simpler and leaves no
credential lying around.

## 2. Tag and archive for Nature

Nature expects a DOI for the archived code, not just a GitHub link:

1. Sign in to [Zenodo](https://zenodo.org) with GitHub and enable the `CAD-scV2E2G` repository.
2. `git tag -a v1.0.0 -m "Nature publication release" && git push --tags`, then create a GitHub release
   from that tag. Zenodo mints a DOI automatically.
3. Put the DOI in `CITATION.cff` and in the Code Availability statement.

Do this after the remaining original scripts land, so the archived version is the complete one.

## 3. Before making it public

Work through the release checklist in [`STATUS.md`](STATUS.md). The two items that matter most: replace
the draft dynamic caQTL script, and add the mosaic analysis code for Figure 4.

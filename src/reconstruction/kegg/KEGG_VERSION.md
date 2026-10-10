kegg118

This is the raven-data release tag that every KEGG-artefact download in
this folder reads from --- getModelFromKEGG's core reference-model/table
bundle, and getKEGGModelForOrganism's HMM libraries. It is the single
place this version is recorded; every function that downloads a KEGG
artefact reads it from here via keggDataVersion.m, instead of
hardcoding a "kegg###" string of its own. Bump this file, not the
calling code, when raven-data publishes a new KEGG release.

The version is read as the file's first line; everything below this
point is only for a human reader.

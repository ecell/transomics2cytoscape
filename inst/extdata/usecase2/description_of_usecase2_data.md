## s2t_name_updated.tsv

The trans-omic interaction file of connection information from the top layer Insulin signal to the 2nd transcription factor layer.
The 3rd column is the KEGG ID of the signal molecule, the 6th column is the TRANSFAC ID of the transcription factor motif
(the last column is the name of the transcription factor),
and the 7th column is the trans-omic interaction type.

## t2e_name_updated.tsv

The trans-omic interaction file of connection information from the 2nd transcription factor layer to the 3rd enzyme layer.
The 3rd column is the TRANSFAC ID of the transcription factor motif
(the last column is the name of the transcription factor),
the 6th column is the KEGG ID of enzyme gene,
and the 7th column is the trans-omic interaction type.

## e2r_updated.tsv
The trans-omic interaction file of connection information from the 3rd enzyme layer to the 4th reaction layer.
The 3rd column is the KEGG ID of enzyme gene, the 6th column is the EC number to be catalyzed,
and the 7th column is the trans-omic interaction type.

## m2r_updated.tsv
The trans-omic interaction file of connection information from the 5th metabolite layer to the 4th reaction layer.
The 3rd column is the KEGG ID of metabolite,
the 6th column is the EC number regulated by metabolite-protein interaction,
and the 7th column is the trans-omic interaction type.

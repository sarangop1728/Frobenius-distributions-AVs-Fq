# Frobenius distributions of abelian varieties over finite fields
Github repository to acompany the paper [Frobenius distributions of low dimensional abelian varieties over finite fields](https://arxiv.org/abs/2306.02237), written by Santiago Arango-Piñeros, Deewang Bhamidipati, and Soumya Sankar.
This code is made available so that the reader can replicate our histograms.

To generate some histograms follow these steps: 
1. Clone this repo.
2. Open the terminal and go to the directory that you just cloned.
3. Fire up sage.
4. Type: `load('sage generate_data.sage')`
5. Type: `histogram_from_label('label')`, where `label` is the LMFDB label of the isogeny class you want.
# Frobenius distributions of abelian varieties over finite fields

Github repository to acompany the paper [Frobenius distributions of low dimensional abelian varieties over finite fields](https://arxiv.org/abs/2306.02237), written by Santiago Arango-Piñeros, Deewang Bhamidipati, and Soumya Sankar.
This code is made available so that the reader can replicate our histograms. You will need [SageMath](https://www.sagemath.org/) installed in your computer.

To generate some histograms you can follow these steps. 
1. Clone this repository.
2. Open the terminal and go to the directory that you just cloned.
3. Fire up SageMath.
4. Type: `load('sage generate_histograms.sage')`
5. Type: `histogram_from_label(label)`, where `label` (a string) is the [LMFDB label](https://www.lmfdb.org/Variety/Abelian/Fq/Labels) of the isogeny class you want.
    - This will generate a .pdf histogram with $16^5$ data points, following the conventions of the paper (see the paragraph before Figure 1). This should take aproximately 30 seconds for an isogeny class of dimension < 4.
    - Alternatively, you can type `histogram_from_label(label, n=6, extension='.png')` if you want $16^6$ data points in .png format. This should take about 10 minutes for an isogeny class of dimension < 4. 
6. To generate all the histograms from the paper, type: `paper_histograms(paper_labels, exponents=[4,5])`.
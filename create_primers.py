"""Script for creating gene assembly primers.
    
Primers containing an ambiguous codon at their center are tiled along a gene.

Jesse Bloom, 2013. 

Edited by Adam Dingens Nov 2015 to generate primers of differing lengths to all have a Tm of ~60C. Notes on edits below. 
This script first makes an ORIGINAL primer of specified length (default 37 bps). 
If the ORIGINAL primer has a Tm of greater than MaxPrimerTm, then nucleotides are trimmed off one by one (first 5', then 3', then 5' etc) until the Tm is less than MaxPrimerTm. Note that this could be over a degree less than the MaxPrimerTm. 
If the ORIGINAL primer has a Tm less than MinPrimerTm, then nucelotides are added one by one (first 3', then 5', then 3' etc) until the Tm is over MinPrimerTm. Note that this could be over a degree more than the MinPrimerTm
If the ORIGINAL primer has a Tm of less than MaxPrimerTm but greater than MinPrimerTm, it is not altered. 
The primers are constrained to be between MinPrimerlength and MaxPrimerLength bps long. The Tm of some MaxPrimerLength primers may not be > MinPrimerTemp, and the Tm of some MinPrimerLength primers may not be < MaxPrimerTm.

For command line arguments, run::

    python create_primers.py -h

The  Tm_NN command of the MeltingTemp Module of Biopython (http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html) is used to calculate Tm of primers. 
This calculation is based on nearest neighbor thermodynamics. nucelotides labeled N are given average values in the Tm calculation. 
It is possible to vary salt concentration and other addatives if needed.

Edited by Kate Crawford January 2021 to include options for `NNG` and `NNC`
ambiguous codons to simulate `NNS` mutagenesis as IDT oPools only allow for
`N` or `K` ambiguous nucelotides.
"""


import os
import sys
import math
import re
import argparse
import streamlit as st
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from io import StringIO



def Parser():
    """Returns command line parser."""
    parser = argparse.ArgumentParser(
            description='Script by Adam Dingens and Jesse Bloom to design codon tiling primers with specific melting temperatures.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            )

    parser.add_argument('sequencefile', help="the name of a file giving the sequence for which we are designing the primers. This file should only contain the sequence, and should not have any headers or other content. For the sequence, make the 5' and 3' ends that you do not want to mutate in lower case. Make the portion of the coding sequence that you want to tile with ambiguous codons in upper case. Typically, for example, the first site you mutate would be the initial start codon, so the first *ATG* would be the first upper case letter. You must have at least *(startprimerlength - 3) / 2* nucleotides in lower case at each end of the upper case sequence that you are mutating. This is because at least this much flanking sequence is needed to design primers of the indicated length; more sequence may be required if the primer at either end is extended beyond the startprimerlength.")
    parser.add_argument('primerprefix', help="prefix name to be added to each primer")
    parser.add_argument('firstcodon', type=int, help='number to assign to first codon in infile to mutagenize, used for primer naming')
    parser.add_argument('outfile', help='name of primer output file')
    parser.add_argument('--startprimerlength', type=int, help='starting primer length', default=37)
    parser.add_argument('--maxprimertm', type=float, help="Upper temperature limit for primers.", default=61)
    parser.add_argument('--minprimertm', type=float, help="Lower temperature limit for primers.", default=60)
    parser.add_argument('--minlength', type=int, help='Minimum primer length', default=25)
    parser.add_argument('--maxlength', type=int, help='Maximum primer length', default=51)
    parser.add_argument('--ambiguous_codon', choices={'NNN', 'NNS', 'NNK', 'NNC', 'NNG'},
                        default='NNN', help='What ambiguous codon to use'),
    parser.add_argument('--output', choices={'plates', 'opools'}, default='plates',
                        help='Format the final csv for ordering in plates or oligo pools.')

    return parser


def ReverseComplement(seq):
    """Returns reverse complement of sequence. Preserves upper/lower case."""
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'N':'N', 'n':'n', 'S':'S', 's':'s', 'K':'M', 'k':'m'}
    rc = [d[nt] for nt in seq]
    rc.reverse()
    return ''.join(rc)


def CreateMutForOligosVarLength(seq, primerlength, prefix, firstcodon, maxprimertm, minprimertm, maxlength, minlength, ambiguous_codon):
    """Creates oligos to tile a gene and introduce ambiguous codon at each site.

    *seq* : sequence of the gene. The gene itself should be upper case. The
    flanking regions and start / stop codons should be lower case. 
    All upper case codons are randomized. The length of the lower
    case sequences at each end must be >= (primerlength - 3) / 2.0

    *primerlength* : length of primers. Must be an odd number, so that equal length
    flanking on each side.

    *prefix* : string prefix attached to primer names.

    *firstcodon* : number assigned to first codon in primer name.

    *ambiguous_codon* : ambiguous codon to use.

    Tiles primers across the gene in the forward direction. The primers
    are all of length primerlength with ambiguous codon at the middle codon.
    Note that only upper case letters are randomized.
    Primers are named as follows:

    "%s-for-mut%d" % (prefix, i) -> 5' tiling primers, where i = 2, 3, ...
    In other words, the indices cover codons 2 and up.

    Returns a list of all these primers as *(name, sequence)* 2-tuples.
    """
    n = len(seq)
    assert primerlength % 2 == 1, "primer length not odd"
    initial_flanklength = (primerlength - 3) // 2
    upperseq = ''.join([nt for nt in seq if nt.istitle()])
    assert upperseq in seq, "upper case nucleotides not substring"
    assert len(upperseq) % 3 == 0, "length of upper case not multiple of 3"
    startupper = seq.index(upperseq)
    if startupper < initial_flanklength:
        raise ValueError("not enough 5' lower case flanking nucleotides")
    if n - len(upperseq) - startupper < initial_flanklength:
        raise ValueError("not enough 3' lower case flanking nucleotides")
    ncodons = len(upperseq) // 3
    primers = []
    for icodon in range(ncodons):
        i = startupper + icodon * 3
        primer = "%s%s%s" % (seq[i - initial_flanklength : i], ambiguous_codon, seq[i + 3 : i + 3 + initial_flanklength]) 
        name = "%s-for-mut%d" % (prefix, firstcodon + icodon)
        primerseq = Seq(primer)
        
        TmNN = ('%0.2f' % mt.Tm_NN(primerseq, strict=False)) 
        add_3 = True
        minus_5 = True
        flank5 = flank3 = initial_flanklength
        if float(TmNN) > float(maxprimertm):
            while float(TmNN) > float(maxprimertm) and len(primer) > minlength:
                if minus_5:
                    flank5 -= 1
                    primer = "%s%s%s" % (seq[i - (flank5) : i], ambiguous_codon, seq[i + 3 : i + 3 + flank3])
                    minus_5 = False
                else:
                    flank3 -= 1
                    primer =  "%s%s%s" % (seq[i - (flank5) : i], ambiguous_codon, seq[i + 3 : i + 3 + flank3])

                    minus_5 = True
                if startupper < flank5:
                    raise ValueError("not enough 5' lower case flanking nucleotides")
                if n - len(upperseq) - startupper < flank3:
                    raise ValueError("not enough 3' lower case flanking nucleotides") #not sure if this is correct!

                
                primerseq = Seq(primer)
                TmNN = ('%0.2f' % mt.Tm_NN(primerseq, strict=False)) 
                primerlength= len(primer)
        else: 
            if float(TmNN) < float(minprimertm):
                while float(TmNN) < float(minprimertm) and len(primer) < maxlength:
                    if add_3:
                        flank3 += 1
                        primer = "%s%s%s" % (seq[i - (flank5) : i], ambiguous_codon, seq[i + 3 : i + 3 + flank3])
                        add_3 = False
                    else:
                        flank5 +=1
                        primer = "%s%s%s" % (seq[i - (flank5) : i], ambiguous_codon, seq[i + 3 : i + 3 + flank3])
                        add_3 = True    
                    primerseq = Seq(primer)
                    TmNN = ('%0.2f' % mt.Tm_NN(primerseq, strict=False)) 
                    primerlength= len(primer)
                    if startupper < flank5:
                        raise ValueError("not enough 5' lower case flanking nucleotides")
                    if n - len(upperseq) - startupper < flank3:
                        raise ValueError("not enough 3' lower case flanking nucleotides") 

            else:
                pass
        primers.append((name, primer))
    #print primers 
    return primers



def main():
    st.title("Primer Design Tool")
    st.write("Use this tool to design mutation primers based on your inputs.")

    # Inputs with defaults
    primerlength = st.number_input("Primer Length (Must be > 3 and Odd)", min_value=1, value=41, step=1)
    sequence = st.text_area("Enter DNA Sequence (no spaces or line breaks)", "attacgtgtttacgaagcaaaagctaaaaccaggagctatttaATGGCAACAGTTAACCAGCTGGTACGCAAACCACGTGCTCGCAAAGTTGCGAAAAGCAACGTGCCTGCGCTGGAAGCATGCCCGCAAAAACGTGGCGTATGTACTCGTGTATATACTACCACTCCTAAAAAACCGAACTCTGCGCTGCGTAAAGTATGCCGTGTTCGTCTGACTAACGGTTTCGAAGTGACTTCCTACATCGGTGGTGAAGGTCACAACCTGCAGGAGCACTCCGTGATCCTGATCCGTGGCGGTCGTGTTAAAGACCTCCCGGGTGTTCGTTACCACACCGTACGTGGTGCGCTTGACTGCTCCGGCGTTAAAGACCGTAAGCAGGCTCGTTCCAAGTATGGCGTGAAGCGTCCTAAGGCTTAAtggttctccgttaagtaaggccaaacgttttaacttaaatgtcaaact")
    outfile = st.text_input("Output File Name", "output.txt")
    primerprefix = st.text_input("Primer Prefix", "rpsl1")
    firstcodon = st.number_input("First Codon Number", min_value=1, value=42, step=1)
    maxprimertm = st.number_input("Maximum Primer Tm", min_value=0.0, value=60.0, step=0.1)
    minprimertm = st.number_input("Minimum Primer Tm", min_value=0.0, value=50.0, step=0.1)
    maxlength = st.number_input("Maximum Primer Length", min_value=1, value=60, step=1)
    minlength = st.number_input("Minimum Primer Length", min_value=1, value=40, step=1)
    ambiguous_codon = st.text_input("Ambiguous Codon", "NNK")
    output = st.selectbox("Output Format", ["plates", "opools"], index=0)

    # Validations
    if st.button("Validate and Run"):
        if primerlength <= 3 or primerlength % 2 == 0:
            st.error("Primer length must be greater than 3 and odd.")
            return

        if not sequence.strip():
            st.error("Please enter a valid DNA sequence.")
            return

        sequence = sequence.replace(' ', '').replace('\n', '')
        st.write(f"Read a sequence of length {len(sequence)}.")

        # Primer design logic
        st.write(f"Using primer prefix: {primerprefix}, and first codon number: {firstcodon}.")
        mutforprimers = CreateMutForOligosVarLength(sequence, primerlength, primerprefix, firstcodon, maxprimertm, minprimertm, maxlength, minlength, ambiguous_codon)
        mutrevprimers = [(name.replace('for', 'rev'), ReverseComplement(seq)) for (name, seq) in mutforprimers]
        primers = mutforprimers + mutrevprimers

        st.write(f"Designed {len(mutforprimers)} forward primers and {len(mutrevprimers)} reverse primers.")
        st.write(f"Total primers: {len(primers)}.")

        output_buffer = StringIO()
        # Writing to output buffer
        try:
            if output == 'plates':
                iplate = 1
                for primers_set in [mutforprimers, mutrevprimers]:
                    output_buffer.write(f"\nPlate {iplate}\n")
                    n_in_plate = 0
                    for (name, primer) in primers_set:
                        output_buffer.write(f"{name}, {primer}\n")
                        n_in_plate += 1
                        if n_in_plate == 96:
                            n_in_plate = 0
                            iplate += 1
                            output_buffer.write(f"\nPlate {iplate}\n")
                    if n_in_plate:
                        iplate += 1
            elif output == 'opools':
                output_buffer.write("Pool name,Primer name,Ambiguous codon,Sequence\n")
                for primers_set in [mutforprimers, mutrevprimers]:
                    pool = f"{primerprefix}_{'ForPool' if primers_set == mutforprimers else 'RevPool'}"
                    for (name, primer) in primers_set:
                        output_buffer.write(f"{pool},{name},{ambiguous_codon},{primer}\n")

            st.success("Primer design completed! You can download the output file below.")
            output_buffer.seek(0)

            # Provide a download button for the generated file
            st.download_button(
                label="Download Primer File",
                data=output_buffer.getvalue(),
                file_name="output.txt",
                mime="text/plain",
            )
        except Exception as e:
            st.error(f"Error generating output: {e}")



if __name__ == "__main__":
    main()
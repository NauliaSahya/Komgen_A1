import streamlit as st
import pandas as pd
import random

Codon_DNA = {
    "Stop Codon": {
        "Codon": ['TAA', 'TAG', 'TGA'],
        "Single_Letter": "Stop",
        "3Letter": "Termination"
    },
    "Isoleucine": {
        "Codon": ['ATT', 'ATC', 'ATA'],
        "Single_Letter": "I",
        "3Letter": "Ile"
    },
    "Leucine": {
        "Codon": ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
        "Single_Letter": "L",
        "3Letter": "Leu"
    },
    "Valine": {
        "Codon": ['GTT', 'GTC', 'GTA', 'GTG'],
        "Single_Letter": "V",
        "3Letter": "Val"
    },
    "Phenylalanine": {
        "Codon": ['TTT', 'TTC'],
        "Single_Letter": "F",
        "3Letter": "Phe"
    },
    "Methionine": {
        "Codon": ['ATG'],
        "Single_Letter": "M",
        "3Letter": "Met (Start)"
    },
    "Cysteine": {
        "Codon": ['TGT', 'TGC'],
        "Single_Letter": "C",
        "3Letter": "Cys"
    },
    "Alanine": {
        "Codon": ['GCT', 'GCC', 'GCA', 'GCG'],
        "Single_Letter": "A",
        "3Letter": "Ala"
    },
    "Glycine": {
        "Codon": ['GGT', 'GGC', 'GGA', 'GGG'],
        "Single_Letter": "G",
        "3Letter": "Gly"
    },
    "Proline": {
        "Codon": ['CCT', 'CCC', 'CCA', 'CCG'],
        "Single_Letter": "P",
        "3Letter": "Pro"
    },
    "Threonine": {
        "Codon": ['ACT', 'ACC', 'ACA', 'ACG'],
        "Single_Letter": "T",
        "3Letter": "Thr"
    },
    "Serine": {
        "Codon": ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        "Single_Letter": "S",
        "3Letter": "Ser"
    },
    "Tyrosine": {
        "Codon": ['TAT', 'TAC'],
        "Single_Letter": "Y",
        "3Letter": "Tyr"
    },
    "Tryptophan": {
        "Codon": ['TGG'],
        "Single_Letter": "V",
        "3Letter": "Trp"
    },
    "Glutamine": {
        "Codon": ['CAA', 'CAG'],
        "Single_Letter": "Q",
        "3Letter": "Gln"
    },
    "Asparagine": {
        "Codon": ['AAT', 'AAC'],
        "Single_Letter": "N",
        "3Letter": "Asn"
    },
    "Histidine": {
        "Codon": ['CAT', 'CAC'],
        "Single_Letter": "H",
        "3Letter": "His"
    },
    "Glutamic acid": {
        "Codon": ['GAA', 'GAG'],
        "Single_Letter": "E",
        "3Letter": "Glu"
    },
    "Aspartic acid": {
        "Codon": ['GAT', 'GAC'],
        "Single_Letter": "D",
        "3Letter": "Asp"
    },
    "Lysine": {
        "Codon": ['AAA', 'AAG'],
        "Single_Letter": "K",
        "3Letter": "Lys"
    },
    "Arginine": {
        "Codon": ['CGT', 'CGC', 'CGA', 'AGA', 'AGG'],
        "Single_Letter": "R",
        "3Letter": "Arg"
    }
}

def generate_random_sequence(sequence_type, nitrogen_base, num_sequences):
    random_sequence = []
    for i in range(0, num_sequences):
        number = random.randint(0, len(nitrogen_base) - 1)
        if sequence_type == "RNA" and number == 3:
            number = 4
        elif sequence_type == "DNA" and number == 4:
            number = 3
        random_sequence.append(nitrogen_base[number])
    return random_sequence
    
def split_sequence_into_codons(random_sequence):
    split_sequence = []
    current_string = ""
    for base in random_sequence:
        current_string += base
        if len(current_string) == 3:
            split_sequence.append(current_string)
            current_string = ""
    return split_sequence

def map_codons_to_amino_acids(split_sequence, Codon_DNA):
    for base in split_sequence:
        for amino_acid, data in Codon_DNA.items():
            if base in data["Codon"]:
                st.write(f"Codon {base}: {amino_acid} ({data['Single_Letter']})")

def main():
    purines = ["A", "G"]
    pyrimidines = ["C", "T", "U"]
    nitrogen_base = purines + pyrimidines
    st.sidebar.title("🔬 Genetic Code Explorer")

    selected_option = st.sidebar.selectbox("Choose an option", ["Amino Acid Table 🧬", "Analyze DNA/RNA Sequence 🧪"])

    if selected_option == "Amino Acid Table 🧬":
        st.title("🧬 Amino Acid & Codon Table")
        with st.expander("🔍 Nitrogen Base Classification"):
            st.markdown(f"- **Purines**: {', '.join(purines)}")
            st.markdown(f"- **Pyrimidines**: {', '.join(pyrimidines)}")
            st.markdown("- **DNA**: Adenine (A), Guanine (G), Cytosine (C), Thymine (T)")
            st.markdown("- **RNA**: Adenine (A), Guanine (G), Cytosine (C), Uracil (U)")
        rb = st.radio("Select Table Type:", ["DNA Codon", "RNA Codon"], horizontal=True)
        if rb == "DNA Codon":
            st.subheader("DNA Codon Table with Amino Acids")
        elif rb == "RNA Codon":
            st.subheader("RNA Codon Table with Amino Acids")
        df = pd.DataFrame([(amino_acid, data['Single_Letter'], data['3Letter'], ', '.join(codon.replace("T", "U") if rb == "RNA Codon" else codon for codon in data['Codon'])) for amino_acid, data in Codon_DNA.items()],
                            columns=['Amino Acid', 'Single Letter Code', '3-Letter Code', 'Codons'])
        st.table(df)

    elif selected_option == "Analyze DNA/RNA Sequence 🧪":
        st.title("🧪 DNA/RNA Sequence Analyzer")

        choose = st.selectbox("Select Sequence Type:", ("DNA", "RNA"))
        input_option = st.radio("Select input method:", ("Random", "Manual Input"))

        if input_option == "Random":
            num_sequences = st.number_input("Enter number of bases (must be multiple of 3):", 
                                            min_value=3, step=3)
            random_sequence = generate_random_sequence(choose, nitrogen_base, num_sequences)
        else:
            user_input = st.text_input("Enter your DNA/RNA sequence:")
            # Input Validation
            if user_input:
                user_input = user_input.upper().strip()  # Pastikan uppercase dan tidak ada spasi
                if choose == "DNA" and "U" in user_input:
                    st.error("Invalid DNA sequence! DNA should not contain Uracil (U).")
                elif choose == "RNA" and "T" in user_input:
                    st.error("Invalid RNA sequence! RNA should not contain Thymine (T).")
                elif len(user_input) % 3 != 0:
                    st.error("Sequence length must be a multiple of 3.")
                else:
                    random_sequence = list(user_input)
        
        if st.button("Run"):
            if not random_sequence:
                st.warning("Please provide a valid sequence first!")
            else:
                st.subheader("Generated/Provided Sequence")
                st.write("".join(random_sequence))

                # DNA to RNA 
                if choose == "RNA":
                    Codon_Table = {amino_acid: {"Codon": [codon.replace("T", "U") for codon in data["Codon"]], 
                                                "Single_Letter": data["Single_Letter"]}
                                for amino_acid, data in Codon_DNA.items()}
                else:
                    Codon_Table = Codon_DNA

                # Split sequence
                split_sequence = split_sequence_into_codons(random_sequence)
                st.subheader("Split Sequence into Codons")
                st.write(split_sequence)

                # Mapping Codons to Amino Acids
                st.subheader("Codon to Amino Acid Mapping")
                map_codons_to_amino_acids(split_sequence, Codon_Table)
        
    

    
if __name__ == "__main__":
    main()

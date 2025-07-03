import streamlit as st
import pandas as pd
import random
import plotly.express as px
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np
from io import StringIO
# st.set_page_config(layout="wide")


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
        "Codon": ['CGT', 'CGC', 'CGG', 'CGA', 'AGA', 'AGG'],
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

def count_each_nitrogen_base(sequence):
    count = {"A": 0, "G": 0, "C": 0, "T": 0, "Error": 0}
    for base in sequence:
        if base in count:
            count[base] += 1
        else:
            count["Error"] += 1
    total_bases = sum(count.values())
    count["Total"] = total_bases

    percent = {}
    for base in count:
        if base != "Total": 
            percent[base] = round((count[base] / total_bases) * 100, 2)
    percent["Total"] = round(sum(percent.values()), 1)

    df_seq = pd.DataFrame({
                "Base": ["A", "G", "C", "T", "Error", "Total"], 
                "Count": [count["A"], count["G"], count["C"], count["T"], count["Error"], count["Total"]],
                "Percentage (%)": [percent["A"], percent["G"], percent["C"], percent["T"], percent["Error"], percent["Total"]]
            })

    return count, percent, df_seq

def count_dimer(window, sequence):
    windows = [sequence[i:i+window] for i in range(0, len(sequence), window)]
    
    freq = {
        "A": [], "C": [], "G": [], "T": [],
        "SUM_ACGT": [], "GC": [], "TA": [], "SUM_GC_TA": []
    }
    
    # Hitung frekuensi A,T,G,C, GC content, TA dimer, dan total freq
    for w in windows:
        length = len(w)
        count_A = w.count('A')
        count_C = w.count('C')
        count_G = w.count('G')
        count_T = w.count('T')
        
        freq["A"].append(count_A / length if length > 0 else 0)
        freq["C"].append(count_C / length if length > 0 else 0)
        freq["G"].append(count_G / length if length > 0 else 0)
        freq["T"].append(count_T / length if length > 0 else 0)
        
        # GC content dan TA dimer frekuensi
        gc_content = (count_G + count_C) / length if length > 0 else 0
        
        # Hitung frekuensi dimer TA
        total_dimers = length - 1 if length > 1 else 0
        ta_count = 0
        for i in range(total_dimers):
            if w[i:i+2] == "TA":
                ta_count += 1
        ta_freq = ta_count / total_dimers if total_dimers > 0 else 0
        
        freq["GC"].append(gc_content)
        freq["TA"].append(ta_freq)
        
        freq["SUM_ACGT"].append(freq["A"][-1] + freq["C"][-1] + freq["G"][-1] + freq["T"][-1])
        freq["SUM_GC_TA"].append(gc_content + ta_freq)
    
    return windows, freq


def count_all_dimers_matrix(window, sequence):
    bases = "ATGC"  # urutan sama seperti diminta (A T G C)
    windows = [sequence[i:i+window] for i in range(0, len(sequence), window)]
    
    # Buat dict untuk simpan matriks dimer per window
    dimer_freq_matrices = []
    
    for w in windows:
        length = len(w)
        total_dimers = length - 1 if length > 1 else 0
        
        # Inisialisasi matriks 4x4 dengan 0
        matrix = pd.DataFrame(0, index=list(bases), columns=list(bases), dtype=float)
        
        # Hitung count tiap dimer
        for i in range(total_dimers):
            dimer = w[i:i+2]
            if len(dimer) == 2 and dimer[0] in bases and dimer[1] in bases:
                matrix.at[dimer[0], dimer[1]] += 1
        
        # Konversi count ke frekuensi
        if total_dimers > 0:
            matrix = matrix / total_dimers
        
        dimer_freq_matrices.append(matrix)
    
    return windows, dimer_freq_matrices

def count_basefreq_slide1(win, seq):
    windows = [seq[i:i+win] for i in range(len(seq) - win + 1)]
    freq = {"A": [], "C": [], "G": [], "T": [], "SUM_ACGT": [], "GC": [], "TA": [], "SUM_GC_TA": []}

    for w in windows:
        l = len(w)
        a, c, g, t = w.count('A'), w.count('C'), w.count('G'), w.count('T')
        freq["A"].append(a / l)
        freq["C"].append(c / l)
        freq["G"].append(g / l)
        freq["T"].append(t / l)

        gc = (g + c) / l
        ta = sum(1 for i in range(l - 1) if w[i:i+2] == "TA") / (l - 1) if l > 1 else 0
        freq["GC"].append(gc)
        freq["TA"].append(ta)
        freq["SUM_ACGT"].append((a + c + g + t) / l)
        freq["SUM_GC_TA"].append(gc + ta)

    return windows, freq

def count_dimermat_slide1(win, seq):
    bases = "ATGC"
    windows = [seq[i:i+win] for i in range(len(seq) - win + 1)]
    mats = []

    for w in windows:
        l = len(w)
        total = l - 1
        mat = pd.DataFrame(0, index=list(bases), columns=list(bases), dtype=float)

        for i in range(total):
            d = w[i:i+2]
            if d[0] in bases and d[1] in bases:
                mat.at[d[0], d[1]] += 1

        if total > 0:
            mat /= total

        mats.append(mat)

    return windows, mats

def translate_codon(seq, codon_dict):
    single_letter_seq = ""
    codon_detail = []
    codon = ""
    for base in seq:
        codon += base
        if len(codon) == 3:
            for aa, info in codon_dict.items():
                if codon in info["Codon"]:
                    sl = info["Single_Letter"]
                    single_letter_seq += sl
                    codon_detail.append((codon, aa, sl))
                    break
            codon = ""
    return single_letter_seq, codon_detail

def orf_finder(genome_sequence, Codon_DNA, orf_code, index_orf, min_length, plus_or_min):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    final_dictionary = {}
    index_orf_in_each_code = 1
    
    if plus_or_min == "+": 
        operator_frame = "+"
    elif plus_or_min == "-":
        operator_frame = "-"

    if orf_code == 1:
        i = 0
        type_frame = 1
    elif orf_code == 2:
        i = 1
        type_frame = 2
    elif orf_code == 3:
        i = 2
        type_frame = 3
    
    index_orf += 1
    
    while i < len(genome_sequence):
        if i + 2 < len(genome_sequence) and genome_sequence[i] == 'A' and genome_sequence[i+1] == 'T' and genome_sequence[i+2] == 'G': #find start codon
            j = i + 3 #analyze codon while finding stop codon
            while j < len(genome_sequence) - 2:
                codon = genome_sequence[j] + genome_sequence[j+1] + genome_sequence[j+2]
                if codon in stop_codons: 
                    orf_sequence = genome_sequence[i:j+3] #complete ORF sequence
                    if len(orf_sequence) >= min_length:
                        if operator_frame == "+":
                            amino_acid_seq, codon_details = translate_codon(orf_sequence, Codon_DNA)
                            orf_info = {
                                "Index_ORF" : index_orf,
                                "Frame ORF": f" Frame [{operator_frame}{type_frame}] Number in Each Frame: {index_orf_in_each_code}",
                                "Index Start Codon": i+1,
                                "Index Stop Codon": j+3,
                                "ORF Lenght" :(j+2) - (i+1) +2,
                                "Sequence ORF": orf_sequence,
                                "Sequence Amino Acid": amino_acid_seq,
                                "Codon Detail": codon_details
                            }
                        elif operator_frame == "-":
                            amino_acid_seq, codon_details = translate_codon(orf_sequence, Codon_DNA)
                            orf_info = {
                                "Index_ORF" : index_orf,
                                "Frame ORF": f" Frame [{operator_frame}{type_frame}] Number in Each Frame: {index_orf_in_each_code}",
                                "Index Start Codon": len(genome_sequence) - i,
                                "Index Stop Codon": len(genome_sequence) - j-2,
                                "ORF Lenght" :(j+2) - (i+1) +2,
                                "Sequence ORF": orf_sequence,
                                "Sequence Amino Acid": amino_acid_seq,
                                "Codon Detail": codon_details
                            }
                        final_dictionary[f"Index ORF_{index_orf_in_each_code}"] = orf_info
                        index_orf += 1
                        index_orf_in_each_code += 1
                    i = j + 3 #Move to the position after the stop codon to find the next ORF.
                    break
                j += 3 #continie searching stop codon
            else:
                i += 3 # If no stop codon is found, continue searching at the next position.
        else:
            i += 3 #if no start codon is found, continue searching at the next position.
    
    return final_dictionary

def download_combined_txt(description, genome_sequence, genome_stats, inverse_complement, inv_complement_stats, combined_dicts):
    content = []

    # 1. Deskripsi dan total panjang
    content.append("=== GENOME INFO ===")
    content.append(f"Description: {description}")
    content.append(f"Total Genome Length: {len(genome_sequence)} bp\n")

    # 2. Genome Sequence
    content.append("=== COMPLETE GENOME SEQUENCE ===")
    content.append(genome_sequence + "\n")

    # 3. Statistik Genome
    content.append("=== GENOME STATISTICS ===")
    for base, count in genome_stats[0].items():
        content.append(f"{base}: {count} ({genome_stats[1][base]}%)")
    content.append("")

    # 4. Inverse Complement Genome
    content.append("=== INVERSE COMPLEMENT GENOME SEQUENCE ===")
    content.append(inverse_complement + "\n")

    # 5. Statistik Inverse Complement
    content.append("=== INVERSE COMPLEMENT STATISTICS ===")
    for base, count in inv_complement_stats[0].items():
        content.append(f"{base}: {count} ({inv_complement_stats[1][base]}%)")
    content.append("")

    # 6. ORF List
    total_orf = sum(len(d) for d in combined_dicts)
    content.append(f"=== ORF FOUND: {total_orf} ===\n")
    for idx, dictionary in enumerate(combined_dicts, start=1):
        for orf_key, orf_info in dictionary.items():
            content.append(f"--- ORF {orf_info['Index_ORF']} ---")
            for key, value in orf_info.items():
                if key == "Codon Detail":
                    content.append("Codon Detail:")
                    for codon_idx, (codon, aa_3, aa_1) in enumerate(value, start=1):
                        content.append(f"  codon[{idx}, {codon_idx}]: {codon} {aa_3} ({aa_1})")
                else:
                    content.append(f"{key}: {value}")
            content.append("")

    # Gabung semua jadi satu string
    final_text = "\n".join(content)

    # Tampilkan tombol download
    st.download_button(
        label="Download All as TXT",
        data=final_text,
        file_name="genome_report.txt",
        mime="text/plain"
    )

# === Fungsi Translasi ke Asam Amino ===
def print_codon_translation(obs_raw, path_state, codon_dict, output_txt):
    coding_patterns = [
        ["A_start", "T_start", "G_start"],
        ["Coding1", "Coding2", "Coding3"],
        ["Stop1", "Stop2", "Stop3"]
    ]

    def is_coding_triplet(states):
        return states in coding_patterns

    st.subheader("\n=== CODING CODON TRANSLATION ===")
    output_txt.write("\n=== CODING CODON TRANSLATION ===\n")
    codon_idx = 1
    i = 0
    while i <= len(obs_raw) - 3:
        states = path_state[i:i+3]
        if is_coding_triplet(states):
            codon = obs_raw[i] + obs_raw[i+1] + obs_raw[i+2]
            found = False
            for aa, data in codon_dict.items():
                if codon.upper() in data["Codon"]:
                    line = f"codon[{codon_idx}]: {codon} {aa} ({data['Single_Letter']}) [{i}, {i+1}, {i+2}]"
                    st.write(line)
                    output_txt.write(line + '\n')
                    found = True
                    break
            if not found:
                line = f"codon[{codon_idx}]: {codon} ??? (Unknown) [{i}, {i+1}, {i+2}]"
                st.write(line)
                output_txt.write(line + '\n')
            codon_idx += 1
            i += 3
        else:
            i += 1

    st.subheader("\n=== NON-CODING CODONS ===")
    output_txt.write("\n=== NON-CODING CODONS ===\n")
    noncoding_idx = 1
    j = 0
    while j <= len(obs_raw) - 3:
        states = path_state[j:j+3]
        if not is_coding_triplet(states):
            codon = obs_raw[j] + obs_raw[j+1] + obs_raw[j+2]
            line = f"noncoding[{noncoding_idx}]: {codon} [{j}, {j+1}, {j+2}]"
            st.write(line)
            output_txt.write(line + '\n')
            noncoding_idx +=1
        j += 3


def print_codon_translation8(obs_raw, path_state, codon_dict, output_txt):
    st.subheader("\n=== CODING CODON TRANSLATION ===")
    output_txt.write("\n=== CODING CODON TRANSLATION ===\n")
    codon_idx = 1
    i = 0
    while i < len(path_state):
        # Deteksi awal exon
        if path_state[i] != "non coding" and path_state[i] !="stop1":
            exon_start = i
            exon_seq = []
            while i < len(path_state) and path_state[i] != "non coding" and path_state[i] !="stop1":
                exon_seq.append(obs_raw[i])
                i += 1
            # Proses exon jika panjangnya >= 3
            if len(exon_seq) >= 3:
                trimmed_len = len(exon_seq) - (len(exon_seq) % 3)
                exon_seq = exon_seq[:trimmed_len]
                for j in range(0, len(exon_seq), 3):
                    codon = exon_seq[j:j+3]
                    codon_str = ''.join(codon)
                    found = False
                    for aa, data in codon_dict.items():
                        if codon_str.upper() in data["Codon"]:
                            line =(f"codon[{codon_idx}]: {codon_str} {aa} ({data['Single_Letter']})")
                            st.write(line)
                            output_txt.write(line + '\n')
                            found = True
                            break
                    if not found:
                        line = (f"codon[{codon_idx}]: {codon_str} ??? (Unknown)")
                        st.write(line)
                        output_txt.write(line + '\n')
                    codon_idx += 1
        else:
            i += 1

    st.subheader("\n=== NON-CODING CODONS ===")
    output_txt.write("\n=== NON-CODING CODONS ===\n")
    noncoding_idx = 1
    j = 0
    while j <= len(obs_raw) - 3:
        states_triplet = path_state[j:j+3]
        if all(s == "non coding" for s in states_triplet):
            codon = obs_raw[j:j+3]
            line =(f"noncoding[{noncoding_idx}]: {''.join(codon)} [{j}, {j+1}, {j+2}]")
            st.write(line)
            output_txt.write(line + '\n')
            noncoding_idx += 1
        j += 3

def print_codon_translation12(obs_raw, path_state, codon_dict, output_txt):
    st.subheader("\n=== CODING CODON TRANSLATION ===")
    output_txt.write("\n=== CODING CODON TRANSLATION ===\n")
    codon_idx = 1
    i = 0

    coding_states = ["A_start", "T_start", "G_start", "coding1", "coding2", "coding3", "stop1", "stop2", "stop3"]
    noncoding_states = ["non coding1", "non coding2", "non coding3"]
    # stop_states = ["stop1", "stop2", "stop3"]

    while i < len(path_state):
        if path_state[i] in coding_states:
            exon_seq = []
            exon_start = i
            while i < len(path_state) and path_state[i] in coding_states:
                exon_seq.append(obs_raw[i])
                i += 1

            if len(exon_seq) >= 3:
                trimmed_len = len(exon_seq) - (len(exon_seq) % 3)
                exon_seq = exon_seq[:trimmed_len]
                for j in range(0, len(exon_seq), 3):
                    codon = exon_seq[j:j+3]
                    codon_str = ''.join(codon)
                    found = False
                    for aa, data in codon_dict.items():
                        if codon_str.upper() in data["Codon"]:
                            line = f"codon[{codon_idx}]: {codon_str} {aa} ({data['Single_Letter']})"
                            st.write(line)
                            output_txt.write(line + '\n')
                            found = True
                            break
                    if not found:
                        line = f"codon[{codon_idx}]: {codon_str} ??? (Unknown)"
                        st.write(line)
                        output_txt.write(line + '\n')
                    codon_idx += 1
        else:
            i += 1

    st.subheader("\n=== NON-CODING CODONS ===\n")
    output_txt.write("\n=== NON-CODING CODONS ===\n")
    noncoding_idx = 1
    j = 0
    while j <= len(obs_raw) - 3:
        states_triplet = path_state[j:j+3]
        if all(s in noncoding_states for s in states_triplet):
            codon = obs_raw[j:j+3]
            line = f"noncoding[{noncoding_idx}]: {''.join(codon)} [{j}, {j+1}, {j+2}]"
            st.write(line)
            output_txt.write(line + '\n')
            noncoding_idx += 1
        j += 3

def viterbi(obs_seq, A, B, pi):
    T = len(obs_seq)
    N = A.shape[0]
    delta = np.full((T, N), -np.inf)  # log(0) = -inf
    psi = np.zeros((T, N), dtype=int)

    # Precompute log matrices, ganti 0 dengan nilai kecil agar log tidak error
    logA = np.log(np.where(A == 0, 1e-300, A))
    logB = np.log(np.where(B == 0, 1e-300, B))
    logpi = np.log(np.where(pi == 0, 1e-300, pi))

    # Initialization
    delta[0] = logpi + logB[:, obs_seq[0]]

    # Recursion
    for t in range(1, T):
        for j in range(N):
            prob = delta[t - 1] + logA[:, j] + logB[j, obs_seq[t]]
            delta[t, j] = np.max(prob)
            psi[t, j] = np.argmax(prob)

    # Termination
    path = np.zeros(T, dtype=int)
    path[T - 1] = np.argmax(delta[T - 1])
    for t in range(T - 2, -1, -1):
        path[t] = psi[t + 1, path[t + 1]]

    return path




def main():
    purines = ["A", "G"]
    pyrimidines = ["C", "T", "U"]
    nitrogen_base = purines + pyrimidines
    st.sidebar.title("🔬 Genetic Code Explorer")

    selected_option = st.sidebar.selectbox("Choose an option", ["Run Hidden Markov Model 🧠", "Amino Acid Table 🧬", "Analyze DNA/RNA Sequence 🧪", "Fasta File", "Fasta File with MAV", "Find Open Reading Frame", "ORF Finder per Frame" , "Sequence Alignment🔎"])

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
            num_sequences = st.number_input("Enter number of bases:", 
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
                else:
                    random_sequence = list(user_input)
        
        if st.button("Run"):
            if not random_sequence:
                st.warning("Please provide a valid sequence first!")
            else:
                # 1. Tampilkan urutan asli
                st.subheader("🧬 Generated/Provided Sequence (Original)")
                full_seq = "".join(random_sequence)
                st.write(full_seq)

                # 2. Trimming
                trimmed_length = len(random_sequence) - (len(random_sequence) % 3)
                trimmed_sequence = random_sequence[:trimmed_length]
                if len(random_sequence) % 3 != 0:
                    st.warning(f"Sequence trimmed to {trimmed_length} bases. Last {len(random_sequence) % 3} base(s) ignored.")

                st.subheader("✂️ Trimmed Sequence for Codon Processing")
                st.write("".join(trimmed_sequence))

                # 3. Siapkan codon table
                if choose == "RNA":
                    Codon_Table = {
                        amino_acid: {
                            "Codon": [codon.replace("T", "U") for codon in data["Codon"]],
                            "Single_Letter": data["Single_Letter"]
                        }
                        for amino_acid, data in Codon_DNA.items()
                    }
                else:
                    Codon_Table = Codon_DNA

                # 4. Split menjadi codons
                split_sequence = split_sequence_into_codons(trimmed_sequence)
                st.subheader("🔗 Split Sequence into Codons")
                st.write(split_sequence)

                # 5. Mapping codons to amino acids
                st.subheader("🧬 Codon to Amino Acid Mapping")
                map_codons_to_amino_acids(split_sequence, Codon_Table)


    elif selected_option == "Fasta File":
        st.title("Analyze Fasta File")
        st.subheader("Uploaded FASTA File")
        uploaded_file = st.file_uploader("Upload FASTA File", type=["fasta"])
        genome_sequence = ""
        desc = ""

        if uploaded_file is not None:
            for line in uploaded_file:
                line = line.decode("utf-8").strip()
                if not line.startswith(">"):
                    genome_sequence += line
                else:
                    desc = line
            
            st.subheader("📄 File Information")
            st.write(desc)
            st.write(f"Total Sequence Length: {len(genome_sequence)} bases")

            window = st.number_input("Enter Segment Size", min_value=1, value=10)
            st.subheader("🧪 Sequence Statistic")
            # fungsi count_each_nitrogen_base harus sudah ada sebelumnya
            base_count, base_percent, df = count_each_nitrogen_base(genome_sequence)
            st.write(df)
                
            df_filtered = df[df["Base"].isin(["A", "C", "G", "T"])]
            fig, ax = plt.subplots(figsize=(2, 2))
            ax.pie(df_filtered["Count"], labels=df_filtered["Base"], autopct='%1.2f%%', startangle=90,
                colors=["#FF9999", "#66B3FF", "#99FF99", "#FFCC99"], textprops={'fontsize': 4})
            ax.set_title("Distribusi Nitrogen Base (A, C, G, T)", fontsize=6)
            col1, col2, col3 = st.columns(3)
            with col1:
                st.pyplot(fig)

            # Hitung frekuensi basa dan dimer GC, TA, dll
            windows, freq = count_dimer(window, genome_sequence)
                
            st.subheader("📊 Window Sequence Analysis")
            df_seq = pd.DataFrame({"Window": range(len(windows)), "Sequence": windows})
            df_seq["A Frequency"] = freq["A"]
            df_seq["T Frequency"] = freq["T"]
            df_seq["G Frequency"] = freq["G"]
            df_seq["C Frequency"] = freq["C"]
            df_seq["Total A+C+G+T"] = freq["SUM_ACGT"]
            df_seq["AT Frequency"] = df_seq["A Frequency"] + df_seq["T Frequency"]
            df_seq["GC Frequency"] = df_seq["G Frequency"] + df_seq["C Frequency"]
            df_seq["AT+GC Frequency"] = df_seq["AT Frequency"] + df_seq["GC Frequency"]
                
            st.dataframe(df_seq, use_container_width=True)
                
            fig = go.Figure()

            for base in ["A", "T", "G", "C"]:
                fig.add_trace(go.Scatter(
                    y=freq[base], 
                    mode='lines', 
                    name=f'{base} Frequency'
                ))
            fig.update_layout(
                title='Frequencies of A, T, G, C Content',
                xaxis_title='Window Index',
                yaxis_title='Frequency',
                height=500
            )

            st.plotly_chart(fig, use_container_width=True)

            # Plot GC dan TA frequency (dari freq)
            fig_gc_ta = go.Figure()
            fig_gc_ta.add_trace(go.Scatter(y=freq["GC"], mode='lines', name='GC Content'))
            fig_gc_ta.add_trace(go.Scatter(y=freq["TA"], mode='lines', name='TA Frequency'))
            fig_gc_ta.update_layout(title='GC Content and TA Frequency in Windows',
                                    xaxis_title='Window Index', yaxis_title='Frequency')
            st.plotly_chart(fig_gc_ta, use_container_width=True)

            # Hitung semua frekuensi dimer matriks
            windows2, dimer_freq_matrices = count_all_dimers_matrix(window, genome_sequence)

            st.subheader("🔍 View Dimer Frequency Matrix")

            selected_window = st.number_input(
                "Select window index", 
                min_value=0, 
                max_value=len(dimer_freq_matrices) - 1, 
                value=0, 
                step=1
            )

            matrix = dimer_freq_matrices[selected_window].copy()

            # Tambahkan kolom total per baris
            matrix["Total"] = matrix.sum(axis=1)

            # Tambahkan baris total per kolom
            matrix.loc["Total"] = matrix.sum()

            st.write(f"🪟 Dimer Frequency Matrix - Window {selected_window}")
            st.dataframe(matrix.style.format("{:.4f}"), use_container_width=True)


            # Plot GC+TA dimer
            gc_dimer_freq = []
            ta_dimer_freq = []
            for matrix in dimer_freq_matrices:
                gc = matrix.at['G', 'C'] if 'G' in matrix.index and 'C' in matrix.columns else 0
                ta = matrix.at['T', 'A'] if 'T' in matrix.index and 'A' in matrix.columns else 0
                gc_dimer_freq.append(gc)
                ta_dimer_freq.append(ta)
            fig_dimer_line = go.Figure()
            fig_dimer_line.add_trace(go.Scatter(y=gc_dimer_freq, mode='lines', name='GC Dimer Frequency'))
            fig_dimer_line.add_trace(go.Scatter(y=ta_dimer_freq, mode='lines', name='TA Dimer Frequency'))
            fig_dimer_line.update_layout(title='Dimer Frequencies (GC and TA) per Window',
                                        xaxis_title='Window Index',
                                        yaxis_title='Frequency')
            st.plotly_chart(fig_dimer_line, use_container_width=True)

    elif selected_option == "Fasta File with MAV":
        st.title("Analyze Fasta File")
        st.subheader("Uploaded FASTA File")
        uploaded_file = st.file_uploader("Upload FASTA File", type=["fasta"])
        genome_sequence = ""
        desc = ""

        if uploaded_file is not None:
            for line in uploaded_file:
                line = line.decode("utf-8").strip()
                if not line.startswith(">"):
                    genome_sequence += line
                else:
                    desc = line

            st.subheader("📄 File Information")
            st.write(desc)
            st.write(f"Total Sequence Length: {len(genome_sequence)} bases")

            window = st.number_input("Enter Segment Size", min_value=1, value=10)
            st.subheader("🧪 Sequence Statistic")
            base_count, base_percent, df = count_each_nitrogen_base(genome_sequence)
            st.write(df)

            df_filtered = df[df["Base"].isin(["A", "C", "G", "T"])]
            fig, ax = plt.subplots(figsize=(2, 2))
            ax.pie(df_filtered["Count"], labels=df_filtered["Base"], autopct='%1.2f%%', startangle=90,
                colors=["#FF9999", "#66B3FF", "#99FF99", "#FFCC99"], textprops={'fontsize': 4})
            ax.set_title("Distribusi Nitrogen Base (A, C, G, T)", fontsize=6)
            col1, col2, col3 = st.columns(3)
            with col1:
                st.pyplot(fig)

            # ======== Hitung frekuensi konten dengan step 1 ========
            windows, freq = count_basefreq_slide1(window, genome_sequence)

            st.subheader("📊 Window Sequence Analysis")
            df_seq = pd.DataFrame({"Window": range(len(windows)), "Sequence": windows})
            df_seq["A Frequency"] = freq["A"]
            df_seq["T Frequency"] = freq["T"]
            df_seq["G Frequency"] = freq["G"]
            df_seq["C Frequency"] = freq["C"]
            df_seq["Total A+C+G+T"] = freq["SUM_ACGT"]
            df_seq["AT Frequency"] = df_seq["A Frequency"] + df_seq["T Frequency"]
            df_seq["GC Frequency"] = df_seq["G Frequency"] + df_seq["C Frequency"]
            df_seq["AT+GC Frequency"] = df_seq["AT Frequency"] + df_seq["GC Frequency"]
            st.dataframe(df_seq, use_container_width=True)
            st.write(f"Total Valid Windows: {len(windows)}")

            fig = go.Figure()
            for base in ["A", "T", "G", "C"]:
                fig.add_trace(go.Scatter(
                    y=freq[base],
                    mode='lines',
                    name=f'{base} Frequency'
                ))
            fig.update_layout(
                title='Frequencies of A, T, G, C Content',
                xaxis_title='Window Index',
                yaxis_title='Frequency',
                height=500
            )
            st.plotly_chart(fig, use_container_width=True)

            # Plot GC dan TA frequency (dari freq)
            fig_gc_ta = go.Figure()
            fig_gc_ta.add_trace(go.Scatter(y=freq["GC"], mode='lines', name='GC Content'))
            fig_gc_ta.add_trace(go.Scatter(y=freq["TA"], mode='lines', name='TA Frequency'))
            fig_gc_ta.update_layout(title='GC Content and TA Frequency in Windows',
                                    xaxis_title='Window Index', yaxis_title='Frequency')
            st.plotly_chart(fig_gc_ta, use_container_width=True)

            # ======== Hitung semua frekuensi dimer matrix dengan step 1 ========
            windows2, dimer_freq_matrices = count_dimermat_slide1(window, genome_sequence)

            st.subheader("🔍 View Dimer Frequency Matrix")
            selected_window = st.number_input(
                "Select window index",
                min_value=0,
                max_value=len(dimer_freq_matrices) - 1,
                value=0,
                step=1
            )
            matrix = dimer_freq_matrices[selected_window].copy()
            matrix["Total"] = matrix.sum(axis=1)
            matrix.loc["Total"] = matrix.sum()
            st.write(f"🪟 Dimer Frequency Matrix - Window {selected_window}")
            st.dataframe(matrix.style.format("{:.4f}"), use_container_width=True)

            # ======== Plot GC + TA Dimer per window (line chart) ========
            gc_dimer_freq = []
            ta_dimer_freq = []
            for mat in dimer_freq_matrices:
                gc = mat.at['G', 'C'] if 'G' in mat.index and 'C' in mat.columns else 0
                ta = mat.at['T', 'A'] if 'T' in mat.index and 'A' in mat.columns else 0
                gc_dimer_freq.append(gc)
                ta_dimer_freq.append(ta)

            fig_dimer_line = go.Figure()
            fig_dimer_line.add_trace(go.Scatter(y=gc_dimer_freq, mode='lines', name='GC Dimer Frequency'))
            fig_dimer_line.add_trace(go.Scatter(y=ta_dimer_freq, mode='lines', name='TA Dimer Frequency'))
            fig_dimer_line.update_layout(title='Dimer Frequencies (GC and TA) per Window',
                                        xaxis_title='Window Index',
                                        yaxis_title='Frequency')
            st.plotly_chart(fig_dimer_line, use_container_width=True)

    elif selected_option == "Find Open Reading Frame":
        st.title("🧬Find Gene in the DNA Sequence🧬")
        st.subheader("Uploaded FASTA File")
        #Elemen upload file
        uploaded_file = st.file_uploader("Upload FASTA File", type=["fasta"])
        min_length = st.number_input("Given minimal lenght", min_value=1)
        st.write("k_value is a variabel that define the nitrogen base minimum lenght of the ORF")
        genome_sequence = ""
        desc = ""
        if st.button("Run"):
            for line in uploaded_file:
                line = line.decode("utf-8").strip()
                if not line.startswith(">"):
                    genome_sequence += line
                else:
                    desc = line

            st.subheader("Description of the Genom")
            st.write(desc)
            st.write(f"Total Sequence Length: {len(genome_sequence)} bases")

            st.write("**🧬Complete Genome Sequence:**") #Complete Coding Sequence
            st.write(genome_sequence)

            st.write("**🧪 Sequence Statistic:**")
            base_count, base_percent, df = count_each_nitrogen_base(genome_sequence)
            st.write(df)

            # inverse complement
            invcom_genome_sequence = ""
            for nucleotide in reversed(genome_sequence):
                if nucleotide == "A":
                    invcom_genome_sequence += "T"
                elif nucleotide == "T":
                    invcom_genome_sequence += "A"
                elif nucleotide == "C":
                    invcom_genome_sequence += "G"
                elif nucleotide == "G":
                    invcom_genome_sequence += "C"
            
            st.write("**🧬Inverse Complement Genome Sequence:**")
            st.write(invcom_genome_sequence)

            st.write("**🧪 Sequence Statistic:**")
            invc_count, invc_percent, df = count_each_nitrogen_base(invcom_genome_sequence)
            st.write(df)

            orf_positif_1 = orf_finder(genome_sequence, Codon_DNA, 1, 0, min_length, "+")
            orf_positif_2 = orf_finder(genome_sequence, Codon_DNA, 2, len(orf_positif_1), min_length, "+")
            orf_positif_3 = orf_finder(genome_sequence, Codon_DNA, 3, len(orf_positif_1) + len(orf_positif_2), min_length, "+")
            orf_negatif_1 = orf_finder(invcom_genome_sequence, Codon_DNA, 1, len(orf_positif_1) + len(orf_positif_2)+ len(orf_positif_3), min_length, "-")
            orf_negatif_2 = orf_finder(invcom_genome_sequence, Codon_DNA, 2, len(orf_negatif_1)+len(orf_positif_1) + len(orf_positif_2)+ len(orf_positif_3), min_length, "-")
            orf_negatif_3 = orf_finder(invcom_genome_sequence, Codon_DNA, 3, len(orf_negatif_2) + len(orf_negatif_1)+len(orf_positif_1) + len(orf_positif_2)+ len(orf_positif_3), min_length, "-")
            combined_dicts = [orf_positif_1, orf_positif_2, orf_positif_3, orf_negatif_1, orf_negatif_2, orf_negatif_3]

            st.subheader("Number ORF Found")

            st.write(f"(ORF +1):   {len(orf_positif_1)}")
            st.write(f"(ORF +2):   {len(orf_positif_2)}")
            st.write(f"(ORF +3):   {len(orf_positif_3)}")
            st.write(f"(ORF -1):   {len(orf_negatif_1)}")
            st.write(f"(ORF -2):   {len(orf_negatif_2)}")
            st.write(f"(ORF -3):   {len(orf_negatif_3)}")
            st.write(f"Total ORF Found:   {len(orf_positif_1)+len(orf_positif_2)+len(orf_positif_3)+len(orf_negatif_1)+len(orf_negatif_2)+len(orf_negatif_3)}")

            st.subheader("ORF List:")
            for dictionary in combined_dicts:
                for orf_index, orf_info in enumerate(dictionary.values(), start=1):
                    st.write("--------")
                    for key, val in orf_info.items():
                        if key == "Codon Detail":
                            st.markdown("**Codon Detail:**")
                            for codon_index, (codon, aa_3letter, aa_1letter) in enumerate(val, start=1):
                                st.markdown(f"codon[{orf_index}, {codon_index}]: `{codon}` {aa_3letter} ({aa_1letter})")
                        else:
                            st.markdown(f"**{key}:** {val}")
        st.write("--------")
        # Fungsi untuk mengunduh data sebagai TXT
        # download_combined_dict(combined_dicts, "combined_dict.txt")
        download_combined_txt(
            description=desc,
            genome_sequence=genome_sequence,
            genome_stats=(base_count, base_percent),
            inverse_complement=invcom_genome_sequence,
            inv_complement_stats=(invc_count, invc_percent),
            combined_dicts=combined_dicts
        )

    elif selected_option == "ORF Finder per Frame":
        st.title("🧬ORF Finder per Frame🧬")
        st.subheader("Uploaded FASTA File")

        uploaded_file = st.file_uploader("Upload FASTA File", type=["fasta"])
        min_length = st.number_input("Minimal Length of ORF", min_value=1)
        frame_option = st.radio(
            "Select Frame",
            options=["+1", "+2", "+3", "-1", "-2", "-3"],
            horizontal=True
        )

        genome_sequence = ""
        desc = ""

        if st.button("Run ORF Frame"):
            for line in uploaded_file:
                line = line.decode("utf-8").strip()
                if not line.startswith(">"):
                    genome_sequence += line
                else:
                    desc = line

            st.subheader("📄 File Information")
            st.write(desc)
            st.write(f"Total Sequence Length: {len(genome_sequence)} bases")
            # st.write("**🧬Complete Genome Sequence:**")
            # st.write(genome_sequence)

            # Statistik urutan
            st.write("**🧪 Sequence Statistic:**")
            base_count, base_percent, df = count_each_nitrogen_base(genome_sequence)
            st.dataframe(df)

            # Generate inverse complement jika diperlukan
            invcom_genome_sequence = ""
            for nucleotide in reversed(genome_sequence):
                if nucleotide == "A":
                    invcom_genome_sequence += "T"
                elif nucleotide == "T":
                    invcom_genome_sequence += "A"
                elif nucleotide == "C":
                    invcom_genome_sequence += "G"
                elif nucleotide == "G":
                    invcom_genome_sequence += "C"

            if frame_option.startswith("-"):
                sequence_to_process = invcom_genome_sequence
                strand = "-"
            else:
                sequence_to_process = genome_sequence
                strand = "+"

            frame_number = int(frame_option[-1])
            plus_or_min = frame_option[0]              # "+" atau "-"
            orf_code = int(frame_option[1])            # 1, 2, atau 3
            index_orf = 0                               

            single_orf_dict = orf_finder(sequence_to_process, Codon_DNA, orf_code, index_orf, min_length, plus_or_min)


            st.subheader(f"Number of ORFs Found in Frame {frame_option}")
            st.write(f"Total ORFs: {len(single_orf_dict)}")

            for orf_index, orf_info in enumerate(single_orf_dict.values(), start=1):
                st.write("--------")
                for key, val in orf_info.items():
                    if key == "Codon Detail":
                        st.markdown("**Codon Detail:**")
                        for codon_index, (codon, aa_3letter, aa_1letter) in enumerate(val, start=1):
                            st.markdown(f"codon[{orf_index}, {codon_index}]: `{codon}` {aa_3letter} ({aa_1letter})")
                    else:
                        st.markdown(f"**{key}:** {val}")

            # === Download TXT Output ===
            output_txt = StringIO()
            output_txt.write("ORF Finder per Frame Output\n\n")
            output_txt.write(f"{desc}\n")
            output_txt.write(f"Frame: {frame_option}\n")
            output_txt.write(f"Total Sequence Length: {len(sequence_to_process)}\n")
            output_txt.write(f"Minimum ORF Length: {min_length}\n")
            output_txt.write(f"Total ORFs Found: {len(single_orf_dict)}\n\n")

            for orf_index, orf_info in enumerate(single_orf_dict.values(), start=1):
                output_txt.write("--------\n")
                for key, val in orf_info.items():
                    if key == "Codon Detail":
                        output_txt.write("Codon Detail:\n")
                        for codon_index, (codon, aa_3letter, aa_1letter) in enumerate(val, start=1):
                            output_txt.write(f"codon[{orf_index}, {codon_index}]: {codon} {aa_3letter} ({aa_1letter})\n")
                    else:
                        output_txt.write(f"{key}: {val}\n")

            st.download_button(
                label="📥 Download ORF Frame TXT",
                data=output_txt.getvalue(),
                file_name=f"orf_frame_{frame_option.replace('-', 'neg_').replace('+', 'pos_')}.txt",
                mime="text/plain"
            )




    elif selected_option == "Sequence Alignment🔎":
        class Blosum50:
            def __init__(self):
                self.residues = "ARNDCQEGHILKMFPSTWYV"
                self.residue_scores = [
                    #  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
                    [5,   -2,  -1, -2, -1,  -1,  -1,  0,   -2,  -1, -2,  -1, -1,  -3, -1,  1,   0,   -3,  -2,  0],
                    [-2,  7,   -1, -2, -4,  1,   0,  -3,   0,  -4, -3,   3,  -2,  -3,  -3,  -1,  -1,  -3,  -1,  -3],
                    [-1,  -1,  7,   2,  -2,  0,   0,  0,   1,  -3, -4,   0,  -2,  -4,  -2,  1,   0,  -4,  -2,  -3],
                    [-2,  -2,  2,   8,  -4,  0,   2,  -1,  -1,  -4, -4,  -1,  -4,  -5,  -1,  0,   -1,  -5,  -3,  -4],
                    [-1,  -4,  -2, -4,  13, -3,  -3,  -3,  -3,  -2, -2,  -3,  -2,  -2,  -4,  -1,  -1,  -5,  -3,  -1],
                    [-1,  1,   0,   0,  -3,  7,   2,  -2,   1,  -3, -2,   2,  0,   -4,  -1,  -1,  0,   -1,  -1,  -3],
                    [-1,  0,   0,   2,  -3,  2,   6,  -3,   0,  -4, -3,   1,  -2,  -3,  -1,  -1,  -1,  -3,  -2,  -3],
                    [0,   -3,  0,   -1, -3,  -2,  -3,  8,   -2,  -4, -4,  -2,  -3,  -4,  -2,  0,   -2,  -3,  -3,  -4],
                    [-2,  0,   1,   -1, -3,  1,   0,  -2,  10,  -4, -3,   0,  -3,  -1,  1,   -1,  -2,  -3,  -3,  -4],
                    [-1,  -4,  -3,  -4,  -2,  -3,  -4,  -4,  -4,  5,   2,   -3,  -1,  -1,  0,   -3,  -1,  -3,  -2,  4],
                    [-2,  -3,  -4,  -4,  -2,  -2,  -3,  -4,  -3,  2,   5,   -3,  2,   1,   -3,  -1,  -1,  -2,  -3,  1],
                    [-1,  3,   0,   -1, -3,  2,   1,  -2,   0,  -3, -3,   6,  -2,  -4,  -1,  -2,  -1,  -2,  -1,  -3],
                    [-1,  -2,  -2,  -4,  -2,  0,   -2,  -3,  -1,  2,   3,   -2,  7,   0,   -3,  -1,  -1,  -3,  -2,  -3],
                    [-3,  -3,  -4,  -5,  -2,  -4,  -3,  -4,  -1,  0,   1,   -4,  0,   8,   -4,  -3,  -4,  -1,  1,   -1],
                    [-1,  -3,  -2,  -1,  -4,  -1,  -1,  -2,  -2,  0,   -3,  -1,  -3,  -4,  10,  -3,  -1,  -4,  -3,  -1],
                    [1,   -1,  1,   0,   -1,  0,   -1,  0,   -1,  -3,  -1,  -2,  -1,  -3,  -3,  5,   2,   -3,  -2,  -3],
                    [0,   -1,  0,   -1,  -1,  -1,  -1,  -2,  -2,  -1,  -1,  -1,  -1,  -2,  -1,  2,   5,   -3,  -2,  0],
                    [-3,  -3,  -4,  -5,  -5,  -1,  -3,  -3,  -3,  -3,  -2,  -3,  -1,  1,   -4,  -3,  -3,  15,  -3,  -3],
                    [-2,  -1,  -2,  -3,  -3,  -1,  -2,  -3,  2,   -1,  -1,  -2,  0,   4,   -3,  -2,  -2,  -3,  8,   2],
                    [0,   -3,  -3,  -4,  -1,  -3,  -3,  -4,  -4,  4,   1,   -3,  -2,  -1,  -1,  -3,  -2,  -3,  -1,  5]
                ]

            def get_score(self, res1, res2):
                if res1 not in self.residues or res2 not in self.residues:
                    return -float('inf')  # atau nilai penalty besar
                idx1 = self.residues.index(res1)
                idx2 = self.residues.index(res2)
                return self.residue_scores[idx1][idx2]


        # Contoh penggunaan
        # blosum50 = Blosum50()
        # score = blosum50.get_score("R", "V")
        # st.write(f"Skor substitusi A dan R: {score}")

        #  Global alignment with the Needleman-Wunsch algorithm (simple gap costs)
        def global_alignment(seq1, seq2, blosum50):
            n = len(seq1)
            m = len(seq2)

            # Matrix skor
            score_matrix = [[0] * (m + 1) for _ in range(n + 1)]

            # Gap penalty
            gap_penalty = -8

            # Inisialization
            for i in range(1, n + 1):
                score_matrix[i][0] = i * gap_penalty
            for j in range(1, m + 1):
                score_matrix[0][j] = j * gap_penalty

            # Fill matrix
            for i in range(1, n + 1):
                for j in range(1, m + 1):
                    match = score_matrix[i - 1][j - 1] + blosum50.get_score(seq1[i - 1], seq2[j - 1])
                    delete = score_matrix[i - 1][j] + gap_penalty
                    insert = score_matrix[i][j - 1] + gap_penalty
                    score_matrix[i][j] = max(match, delete, insert)

            # Traceback
            aligned_seq1 = []
            aligned_seq2 = []
            i, j = n, m
            while i > 0 and j > 0:
                current_score = score_matrix[i][j]
                diagonal_score = score_matrix[i - 1][j - 1]
                up_score = score_matrix[i - 1][j]
                left_score = score_matrix[i][j - 1]

                if current_score == diagonal_score + blosum50.get_score(seq1[i - 1], seq2[j - 1]):
                    aligned_seq1.append(seq1[i - 1])
                    aligned_seq2.append(seq2[j - 1])
                    i -= 1
                    j -= 1
                elif current_score == up_score + gap_penalty:
                    aligned_seq1.append(seq1[i - 1])
                    aligned_seq2.append('-')
                    i -= 1
                else:
                    aligned_seq1.append('-')
                    aligned_seq2.append(seq2[j - 1])
                    j -= 1

            while i > 0:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append('-')
                i -= 1
            while j > 0:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j - 1])
                j -= 1

            aligned_seq1 = ''.join(reversed(aligned_seq1))
            aligned_seq2 = ''.join(reversed(aligned_seq2))
            score = score_matrix[n][m]

            return aligned_seq1, aligned_seq2, score, score_matrix
        
        def parse_fasta(file):
            lines = file.read().decode("utf-8").splitlines()
            desc = ""
            seq = ""
            for line in lines:
                if line.startswith(">"):
                    desc = line[1:].strip()
                else:
                    seq += line.strip()
            return desc, seq.upper()
        
        def create_alignment_marker(aligned1, aligned2):
            marker = ""
            for a, b in zip(aligned1, aligned2):
                if a == b:
                    marker += "|"
                elif a == '-' or b == '-':
                    marker += " "
                else:
                    marker += "."
            return marker
        
        # Local alignment with the Smith-Waterman algorithm (simple gap costs)
        def local_alignment(seq1, seq2, blosum, gap_penalty=-8):
            n, m = len(seq1), len(seq2)
            score_matrix = [[0] * (m + 1) for _ in range(n + 1)]
            max_score = 0
            max_pos = (0, 0)

            # Hitung skor matriks
            for i in range(1, n + 1):
                for j in range(1, m + 1):
                    a, b = seq1[i - 1], seq2[j - 1]
                    subs_score = blosum.get_score(a, b)
                    match = score_matrix[i - 1][j - 1] + subs_score
                    delete = score_matrix[i - 1][j] + gap_penalty
                    insert = score_matrix[i][j - 1] + gap_penalty
                    score_matrix[i][j] = max(0, match, delete, insert)

                    if score_matrix[i][j] > max_score:
                        max_score = score_matrix[i][j]
                        max_pos = (i, j)

            # Traceback
            aligned1, aligned2 = "", ""
            i, j = max_pos
            while i > 0 and j > 0 and score_matrix[i][j] != 0:
                current = score_matrix[i][j]
                a, b = seq1[i - 1], seq2[j - 1]
                subs_score = blosum.get_score(a, b)
                if current == score_matrix[i - 1][j - 1] + subs_score:
                    aligned1 = a + aligned1
                    aligned2 = b + aligned2
                    i -= 1
                    j -= 1
                elif current == score_matrix[i - 1][j] + gap_penalty:
                    aligned1 = a + aligned1
                    aligned2 = "-" + aligned2
                    i -= 1
                elif current == score_matrix[i][j - 1] + gap_penalty:
                    aligned1 = "-" + aligned1
                    aligned2 = b + aligned2
                    j -= 1
                else:
                    break

            return aligned1, aligned2, max_score, score_matrix
        
        # Repeated matches (simple gap costs)
        def repeated_local_alignment(seq1, seq2, blosum, gap_penalty=-8, threshold=0):
            results = []
            used1 = [False] * len(seq1)
            used2 = [False] * len(seq2)

            def mask_sequence(seq, used):
                return ''.join('X' if used[i] else seq[i] for i in range(len(seq)))

            while True:
                masked_seq1 = mask_sequence(seq1, used1)
                masked_seq2 = mask_sequence(seq2, used2)
                n, m = len(masked_seq1), len(masked_seq2)

                # Inisialisasi matriks F
                F = [[0] * (m + 1) for _ in range(n + 1)]

                max_score = 0
                max_pos = (0, 0)

                # Hitung F matrix
                for i in range(1, n + 1):
                    for j in range(1, m + 1):
                        c1, c2 = masked_seq1[i - 1], masked_seq2[j - 1]
                        s = blosum.get_score(c1, c2)
                        match = F[i - 1][j - 1] + s
                        delete = F[i - 1][j] + gap_penalty
                        insert = F[i][j - 1] + gap_penalty
                        F[i][j] = max(0, match, delete, insert)
                        if F[i][j] > max_score:
                            max_score = F[i][j]
                            max_pos = (i, j)

                if max_score <= threshold:
                    break

                # Traceback
                aligned1, aligned2 = "", ""
                i, j = max_pos
                while i > 0 and j > 0 and F[i][j] != 0:
                    c1, c2 = masked_seq1[i - 1], masked_seq2[j - 1]
                    s = blosum.get_score(c1, c2)
                    if F[i][j] == F[i - 1][j - 1] + s:
                        aligned1 = c1 + aligned1
                        aligned2 = c2 + aligned2
                        i -= 1
                        j -= 1
                    elif F[i][j] == F[i - 1][j] + gap_penalty:
                        aligned1 = masked_seq1[i - 1] + aligned1
                        aligned2 = '-' + aligned2
                        i -= 1
                    elif F[i][j] == F[i][j - 1] + gap_penalty:
                        aligned1 = '-' + aligned1
                        aligned2 = masked_seq2[j - 1] + aligned2
                        j -= 1
                    else:
                        break

                # Temukan posisi asli untuk masking
                clean1 = aligned1.replace("-", "")
                clean2 = aligned2.replace("-", "")
                idx1 = masked_seq1.find(clean1)
                idx2 = masked_seq2.find(clean2)

                if idx1 == -1 or idx2 == -1:
                    break

                # Tandai bagian yang telah digunakan
                for k in range(len(clean1)):
                    if idx1 + k < len(used1):
                        used1[idx1 + k] = True
                for k in range(len(clean2)):
                    if idx2 + k < len(used2):
                        used2[idx2 + k] = True

                results.append((aligned1, aligned2, max_score, F))

            return results

        
        # Overlap matching (simple gap costs)
        def overlap_alignment(seq1, seq2, blosum, gap_penalty=-8):
            n = len(seq1)
            m = len(seq2)

            # Initialize score matrix and traceback
            score_matrix = np.zeros((n + 1, m + 1), dtype=int)
            traceback = np.empty((n + 1, m + 1), dtype=object)

            # Initial conditions
            for i in range(1, n + 1):
                score_matrix[i][0] = i * gap_penalty
                traceback[i][0] = 'U'  # Up direction
            for j in range(1, m + 1):
                score_matrix[0][j] = 0  # Overlap, so the first row is 0
                traceback[0][j] = 'L'  # Left direction

            # Fill the score matrix using BLOSUM50 matrix for matches
            for i in range(1, n + 1):
                for j in range(1, m + 1):
                    match = score_matrix[i - 1][j - 1] + blosum.get_score(seq1[i - 1], seq2[j - 1]) if seq1[i - 1] in blosum.residues and seq2[j - 1] in blosum.residues else gap_penalty
                    delete = score_matrix[i - 1][j] + gap_penalty
                    insert = score_matrix[i][j - 1] + gap_penalty
                    score_matrix[i][j] = max(match, delete, insert)

                    # Traceback logic
                    if score_matrix[i][j] == match:
                        traceback[i][j] = 'D'  # Diagonal
                    elif score_matrix[i][j] == delete:
                        traceback[i][j] = 'U'  # Up
                    else:
                        traceback[i][j] = 'L'  # Left

            # Find the maximum score on the last row or column
            max_score = float('-inf')
            i_max, j_max = 0, 0
            for j in range(m + 1):
                if score_matrix[n][j] > max_score:
                    max_score = score_matrix[n][j]
                    i_max, j_max = n, j
            for i in range(n + 1):
                if score_matrix[i][m] > max_score:
                    max_score = score_matrix[i][m]
                    i_max, j_max = i, m

            # Traceback to build the aligned sequences
            aligned1, aligned2 = "", ""
            i, j = i_max, j_max
            while i > 0 and j > 0:
                if traceback[i][j] == 'D':
                    aligned1 = seq1[i - 1] + aligned1
                    aligned2 = seq2[j - 1] + aligned2
                    i -= 1
                    j -= 1
                elif traceback[i][j] == 'U':
                    aligned1 = seq1[i - 1] + aligned1
                    aligned2 = "-" + aligned2
                    i -= 1
                else:  # 'L'
                    aligned1 = "-" + aligned1
                    aligned2 = seq2[j - 1] + aligned2
                    j -= 1

            return aligned1, aligned2, max_score, score_matrix

        # Global alignment using the Needleman-Wunsch algorithm (affine gap costs)
        def affine_global_alignment(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_open=-2, gap_extend=-1):
            n, m = len(seq1), len(seq2)

            # Matriks skor: M = match/mismatch, X = gap di seq1, Y = gap di seq2
            M = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
            X = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
            Y = [[float('-inf')] * (m + 1) for _ in range(n + 1)]
            traceback = [[None] * (m + 1) for _ in range(n + 1)]

            M[0][0] = 0
            for i in range(1, n + 1):
                X[i][0] = gap_open + (i - 1) * gap_extend
                M[i][0] = X[i][0]
            for j in range(1, m + 1):
                Y[0][j] = gap_open + (j - 1) * gap_extend
                M[0][j] = Y[0][j]

            for i in range(1, n + 1):
                for j in range(1, m + 1):
                    match = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty

                    X[i][j] = max(M[i - 1][j] + gap_open + gap_extend,
                                X[i - 1][j] + gap_extend)
                    Y[i][j] = max(M[i][j - 1] + gap_open + gap_extend,
                                Y[i][j - 1] + gap_extend)
                    M[i][j] = max(M[i - 1][j - 1] + match,
                                X[i][j],
                                Y[i][j])

                    if M[i][j] == M[i - 1][j - 1] + match:
                        traceback[i][j] = 'D'
                    elif M[i][j] == X[i][j]:
                        traceback[i][j] = 'U'
                    else:
                        traceback[i][j] = 'L'

            aligned1, aligned2 = "", ""
            i, j = n, m
            while i > 0 or j > 0:
                if i > 0 and j > 0 and traceback[i][j] == 'D':
                    aligned1 = seq1[i - 1] + aligned1
                    aligned2 = seq2[j - 1] + aligned2
                    i -= 1
                    j -= 1
                elif i > 0 and traceback[i][j] == 'U':
                    aligned1 = seq1[i - 1] + aligned1
                    aligned2 = "-" + aligned2
                    i -= 1
                elif j > 0 and traceback[i][j] == 'L':
                    aligned1 = "-" + aligned1
                    aligned2 = seq2[j - 1] + aligned2
                    j -= 1
                else:
                    break

            return aligned1, aligned2, M[n][m], M
        
        # Local alignment with the Smith-Waterman algorithm (affine gap costs, smart linear space algorithm)
        def affine_local_alignment(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_open=-2, gap_extend=-1):
            n, m = len(seq1), len(seq2)

            # Matriks skor
            M = [[0] * (m + 1) for _ in range(n + 1)]
            X = [[0] * (m + 1) for _ in range(n + 1)]
            Y = [[0] * (m + 1) for _ in range(n + 1)]
            traceback = [[None] * (m + 1) for _ in range(n + 1)]

            max_score = 0
            max_pos = (0, 0)

            for i in range(1, n + 1):
                for j in range(1, m + 1):
                    match = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty

                    X[i][j] = max(M[i - 1][j] + gap_open + gap_extend, X[i - 1][j] + gap_extend)
                    Y[i][j] = max(M[i][j - 1] + gap_open + gap_extend, Y[i][j - 1] + gap_extend)
                    M[i][j] = max(0, M[i - 1][j - 1] + match, X[i][j], Y[i][j])

                    if M[i][j] == 0:
                        traceback[i][j] = None
                    elif M[i][j] == M[i - 1][j - 1] + match:
                        traceback[i][j] = 'D'
                    elif M[i][j] == X[i][j]:
                        traceback[i][j] = 'U'
                    else:
                        traceback[i][j] = 'L'

                    if M[i][j] > max_score:
                        max_score = M[i][j]
                        max_pos = (i, j)

            aligned1, aligned2 = "", ""
            i, j = max_pos
            while i > 0 and j > 0 and traceback[i][j]:
                if traceback[i][j] == 'D':
                    aligned1 = seq1[i - 1] + aligned1
                    aligned2 = seq2[j - 1] + aligned2
                    i -= 1
                    j -= 1
                elif traceback[i][j] == 'U':
                    aligned1 = seq1[i - 1] + aligned1
                    aligned2 = "-" + aligned2
                    i -= 1
                elif traceback[i][j] == 'L':
                    aligned1 = "-" + aligned1
                    aligned2 = seq2[j - 1] + aligned2
                    j -= 1
                else:
                    break

            return aligned1, aligned2, max_score, M
        
        # Global alignment (simple gap costs, smart linear-space algorithm)
        def nwsmart_alignment(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-1):
            n, m = len(seq1), len(seq2)

            # Hanya dua kolom matriks yang diperlukan pada waktu yang sama
            prev_col = [0] * (m + 1)
            curr_col = [0] * (m + 1)

            # Mengisi baris pertama
            for j in range(m + 1):
                prev_col[j] = j * gap_penalty

            # Mengisi matriks alignment
            for i in range(1, n + 1):
                curr_col[0] = i * gap_penalty
                for j in range(1, m + 1):
                    match = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty
                    curr_col[j] = max(
                        prev_col[j - 1] + match,     # Diagonal
                        prev_col[j] + gap_penalty,    # Up (gap in seq2)
                        curr_col[j - 1] + gap_penalty # Left (gap in seq1)
                    )
                prev_col, curr_col = curr_col, prev_col  # Berganti kolom

            # Matriks skor untuk output
            score_matrix = pd.DataFrame([prev_col], columns=[f"col_{i+1}" for i in range(m + 1)])

            return prev_col[m], score_matrix
        
        # Local alignment with the Smith-Waterman algorithm (simple gap costs, smart linear space algorithm)
        def swsmart_alignment(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-1):
            n, m = len(seq1), len(seq2)

            # Hanya dua kolom matriks yang diperlukan pada waktu yang sama
            prev_col = [0] * (m + 1)
            curr_col = [0] * (m + 1)

            max_score = 0
            max_pos = (0, 0)

            # Mengisi matriks alignment
            for i in range(1, n + 1):
                curr_col[0] = 0  # Alignment lokal selalu dimulai dari 0
                for j in range(1, m + 1):
                    match = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty
                    curr_col[j] = max(
                        0,  # Zero karena ini alignment lokal
                        prev_col[j - 1] + match,     # Diagonal
                        prev_col[j] + gap_penalty,    # Up (gap in seq2)
                        curr_col[j - 1] + gap_penalty # Left (gap in seq1)
                    )

                    # Menyimpan skor tertinggi dan posisinya
                    if curr_col[j] > max_score:
                        max_score = curr_col[j]
                        max_pos = (i, j)

                prev_col, curr_col = curr_col, prev_col  # Berganti kolom

            # Matriks skor untuk output
            score_matrix = pd.DataFrame([prev_col], columns=[f"col_{i+1}" for i in range(m + 1)])

            return max_score, score_matrix





        
        # Streamlit UI
        st.title("🔬 DNA Sequence Alignment")

        input_mode = st.radio("Pilih metode input:", ["Input", "Upload FASTA"])

        # Default sequence
        default_seq1 = "CCATGGTACATTTGGCTAGGTTTTATAGCTGGCTTGATTGCCATAGTAATGGTGACAATTATGCTTTGCTGTATGACCAGTTGCTGTAGCGC"
        default_seq2 = "CCGTATGTTTGTTAGCAAAACAAGTATCTGTAGATGCTATGTCACGAGTGACACCACCATCAATAGCCTTGTATCCTATGATTTCACTTGACC"

        if input_mode == "Default":
            seq1 = default_seq1
            seq2 = default_seq2
            st.success("Menggunakan sekuens bawaan dari soal.")

        # elif input_mode == "Manual Input":
        #     seq1 = st.text_area("Sequence 1", height=100).upper()
        #     seq2 = st.text_area("Sequence 2", height=100).upper()

        elif input_mode == "Upload FASTA":
            fasta_file1 = st.file_uploader("Upload file FASTA untuk Sequence 1", type=["fasta", "fa"])
            fasta_file2 = st.file_uploader("Upload file FASTA untuk Sequence 2", type=["fasta", "fa"])

            seq1, seq2 = "", ""
            if fasta_file1 and fasta_file2:
                desc1, seq1 = parse_fasta(fasta_file1)
                desc2, seq2 = parse_fasta(fasta_file2)
                st.text(f"Deskripsi Seq1: {desc1}")
                st.text(f"Deskripsi Seq2: {desc2}")

        seq1 = st.text_area("Sequence 1", default_seq1, height=150)
        seq2 = st.text_area("Sequence 2", default_seq2, height=150)
        blosum50 = Blosum50()

        st.write("""
        ### Alignment Marker Explanation:
        - **|**: Indicates a match between the aligned characters.
        - **.**: Indicates a mismatch between the aligned characters.
        - **Space**: Indicates a gap (i.e., '-' in one of the sequences).
        """)

        if st.button("Needleman–Wunsch (Global Alignment)"):
            if not seq1 or not seq2:
                st.error("Pastikan kedua sekuens telah tersedia.")
            else:
                aligned1, aligned2, score, score_matrix = global_alignment(seq1, seq2, blosum50)
                marker = create_alignment_marker(aligned1, aligned2)

                columns = [f"{char}_{i+1}" for i, char in enumerate("-" + seq2)]
                index = [f"{char}_{i+1}" for i, char in enumerate("-" + seq1)]

                score_matrix_df = pd.DataFrame(score_matrix, columns=columns, index=index)

                # Output ke Streamlit
                st.subheader("📊 Matriks Skor (Needleman–Wunsch)")
                st.dataframe(score_matrix_df, height=500)

                # Alignment Result
                st.subheader("🧬 Alignment Result")
                st.text(f"Alignment Score: {score}")
                st.code(f"{aligned1}\n{marker}\n{aligned2}", language="text")
                
        if st.button("Smith–Waterman"):
            if not seq1 or not seq2:
                st.error("Pastikan kedua sekuens telah tersedia.")
            else:
                aligned1, aligned2, score, score_matrix = local_alignment(seq1, seq2, blosum50, gap_penalty=-8)
                marker = create_alignment_marker(aligned1, aligned2)

                # Kolom dan indeks
                columns = [f"{char}_{i}" for i, char in enumerate("-" + seq2)]
                index = [f"{char}_{i}" for i, char in enumerate("-" + seq1)]
                df = pd.DataFrame(score_matrix, columns=columns, index=index)

                # Output
                st.subheader("📊 Matriks Skor (Smith–Waterman)")
                st.dataframe(df, height=500)

                st.subheader("🧬 Alignment Result")
                st.text(f"Alignment Score: {score}")
                st.code(f"{aligned1}\n{marker}\n{aligned2}", language="text")

        if st.button("Repeated Matches"):
            if not seq1 or not seq2:
                st.error("Pastikan kedua sekuens telah tersedia.")
            else:
                results = repeated_local_alignment(seq1, seq2, blosum50)
                if not results:
                    st.warning("Tidak ditemukan local alignment yang signifikan.")
                else:
                    for idx, (a1, a2, score, matrix) in enumerate(results):
                        marker = create_alignment_marker(a1, a2)
                        st.subheader(f"🧬 Alignment {idx+1}")
                        st.text(f"Alignment Score: {score}")
                        st.code(f"{a1}\n{marker}\n{a2}", language="text")
                        df = pd.DataFrame(matrix)
                        st.write("Score Matrix:")
                        st.dataframe(df.style.background_gradient(cmap='Blues', axis=None))
        
        if st.button("Overlap Match"):
            if not seq1 or not seq2:
                st.error("Pastikan kedua sekuens telah tersedia.")
            else:
                aligned1, aligned2, score, score_matrix = overlap_alignment(seq1, seq2, blosum50)
                marker = create_alignment_marker(aligned1, aligned2)

                columns = [f"{char}_{i+1}" for i, char in enumerate("-" + seq2)]
                index = [f"{char}_{i+1}" for i, char in enumerate("-" + seq1)]

                st.subheader("📊 Matriks Skor (Overlap Match)")
                st.dataframe(score_matrix, height=500)

                st.subheader("🧬 Alignment Result")
                st.text(f"Alignment Score: {score}")
                st.code(f"{aligned1}\n{marker}\n{aligned2}", language="text")

        if st.button("Affine Global Alignment"):
            if not seq1 or not seq2:
                st.error("Pastikan kedua sekuens telah tersedia.")
            else:
                aligned1, aligned2, score, score_matrix = affine_global_alignment(seq1, seq2)
                marker = create_alignment_marker(aligned1, aligned2)

                columns = [f"{char}_{i+1}" for i, char in enumerate("-" + seq2)]
                index = [f"{char}_{i+1}" for i, char in enumerate("-" + seq1)]

                st.subheader("📊 Matriks Skor (Affine Global Alignment)")
                st.dataframe(score_matrix, height=500)

                st.subheader("🧬 Alignment Result")
                st.text(f"Alignment Score: {score}")
                st.code(f"{aligned1}\n{marker}\n{aligned2}", language="text")
        
        if st.button("Smart Affine Local Alignment"):
            if not seq1 or not seq2:
                st.error("Pastikan kedua sekuens telah tersedia.")
            else:
                aligned1, aligned2, score, score_matrix = affine_local_alignment(seq1, seq2)
                marker = create_alignment_marker(aligned1, aligned2)

                columns = [f"{char}_{i+1}" for i, char in enumerate("-" + seq2)]
                index = [f"{char}_{i+1}" for i, char in enumerate("-" + seq1)]

                st.subheader("📊 Matriks Skor (Affine Local Alignment)")
                st.dataframe(score_matrix, height=500)

                st.subheader("🧬 Alignment Result")
                st.text(f"Alignment Score: {score}")
                st.code(f"{aligned1}\n{marker}\n{aligned2}", language="text")

        # if st.button("Smart Global Alignment (Linear-Space)"):
        #     # seq1 = st.text_input("Input Sequence 1")
        #     # seq2 = st.text_input("Input Sequence 2")

        #     if not seq1 or not seq2:
        #         st.error("Pastikan kedua sekuens telah tersedia.")
        #     else:
        #         score, score_matrix = nwsmart_alignment(seq1, seq2)

        #         st.subheader("📊 Matriks Skor (Smart Global Alignment)")
        #         st.dataframe(score_matrix, height=500)

        #         st.subheader("🧬 Alignment Result")
        #         st.text(f"Alignment Score: {score}")

        # if st.button("Smart Local Alignment (Linear-Space)"):
        #     # seq1 = st.text_input("Input Sequence 1")
        #     # seq2 = st.text_input("Input Sequence 2")

        #     if not seq1 or not seq2:
        #         st.error("Pastikan kedua sekuens telah tersedia.")
        #     else:
        #         score, score_matrix = swsmart_alignment(seq1, seq2)

        #         st.subheader("📊 Matriks Skor (Smart Local Alignment)")
        #         st.dataframe(score_matrix, height=500)

        #         st.subheader("🧬 Alignment Result")
        #         st.text(f"Alignment Score: {score}")

    elif selected_option == "Run Hidden Markov Model 🧠":
        st.title("🧠 Run Hidden Markov Model (HMM) on FASTA Sequence (Prokaryotic)")
        uploaded_file = st.file_uploader("Upload FASTA file (prokaryotic)", type=["fasta"])
        dat_length = st.number_input("Given lenght to observe (prokaryotic)", min_value=1)
        min_coding_length = st.number_input("Minimal Length for Coding Region", min_value=3, value=15, step=3)
        if st.button("Run Prokaryotic"):
            output_txt = StringIO()
            genome_sequence = ""
            for line in uploaded_file:
                line = line.decode("utf-8").strip()
                if not line.startswith(">"):
                    genome_sequence += line
                else:
                    desc = line   
            st.subheader("📄 File Information")
            output_txt.write("File Information\n")
            st.write(desc)
            output_txt.write(desc + '\n')
            st.write(f"Total Sequence Length: {len(genome_sequence)} bases")
            output_txt.write(f"Total Sequence Length: {len(genome_sequence)} bases\n")

            st.write("**🧬Complete Genome Sequence:**") #Complete Coding Sequence
            output_txt.write("\nComplete Genome Sequence:\n")
            st.write(genome_sequence)
            output_txt.write(genome_sequence)

            st.write("**🧪 Sequence Statistic:**")
            output_txt.write("\n\nSequence Statistic:\n")
            base_count, base_percent, df = count_each_nitrogen_base(genome_sequence)
            st.dataframe(df, use_container_width=False)
            output_txt.write(df.to_string(index=False) + '\n\n')

            # === 1. Definisi State dan Observasi ===
            states = [
                "A_start", "T_start", "G_start", "coding",
                "stop1", "stop2", "stop3", "non coding"
            ]
            n_states = len(states)
            observations = ['A', 'T', 'C', 'G']
            dna_map = {b: i for i, b in enumerate(observations)}

            # === 2. Matriks Transisi ===
            A = np.array([
                [0,   1,   0,   0,   0,   0,   0,   0],
                [0,   0,   1,   0,   0,   0,   0,   0],
                [0,   0,   0,   1,   0,   0,   0,   0],
                [0,   0,   0, 0.9, 0.1,  0,   0,   0],
                [0,   0,   0,   0,   0,   1,   0,   0],
                [0,   0,   0,   0,   0,   0,   1,   0],
                [0,   0,   0,   0,   0,   0,   0,   1],
                [0.1, 0,   0,   0,   0,   0,   0, 0.9]
            ])

            # === 3. Matriks Emisi ===
            B = np.array([
                [1, 0, 0, 0],      # A_start
                [0, 1, 0, 0],      # T_start
                [0, 0, 0, 1],      # G_start
                [0.25, 0.25, 0.25, 0.25],  # Coding
                [0, 1, 0, 0],      # Stop1
                [0.67, 0, 0, 0.33],# Stop2
                [0.67, 0, 0, 0.33],# Stop3
                [0.25, 0.25, 0.25, 0.25],  # NC
            ])

            # === 4. Probabilitas Awal ===
            pi = np.zeros(n_states)
            pi[0] = 0.1  # A_start
            pi[7] = 0.9  # NC

            # === 5. Observasi dari user ===

            # Batasi karakter 
            obs_raw = genome_sequence[:dat_length]

            obs_seq = np.array([dna_map[b] for b in obs_raw])
            
            # === 6. Viterbi Algorithm ===
            path_idx = viterbi(obs_seq, A, B, pi)
            path_state = [states[i] for i in path_idx]

            # === 7. Print Sample Output ===
            # === Print Matriks Transisi A ===
            st.subheader("=== Matriks Transisi (A) ===")
            output_txt.write("\nMatriks Transisi (A)\n")
            df_A = pd.DataFrame(A, columns=states, index=states)
            st.dataframe(df_A, use_container_width=False)
            output_txt.write(df_A.to_string(index=False) + '\n\n')

            # === Print Matriks Emisi B ===
            st.subheader("=== Matriks Emisi (B) ===")
            output_txt.write("\nMatriks Emisi (B)\n")
            df_B = pd.DataFrame(B, columns=["A", "T", "C", "G"], index=states)
            st.dataframe(df_B, use_container_width=False)
            output_txt.write(df_B.to_string(index=False) + '\n\n')


            # === Print Probabilitas Awal π ===
            st.subheader("=== Probabilitas Awal (π) ===")
            output_txt.write("\nProbabilitas Awal (π)\n")
            df_pi = pd.DataFrame({
                "State": states,
                "Initial Probability": pi
            })
            st.dataframe(df_pi, use_container_width=False)
            output_txt.write(df_pi.to_string(index=False) + '\n\n')


            # === Print Semua State dari Viterbi ===
            st.subheader("=== Hidden State Result ===")
            output_txt.write("\nHidden State Result\n")
            viterbi_data = []
            for i, (base, state) in enumerate(zip(obs_raw, path_state)):
                viterbi_data.append({
                    "Position": i,
                    "Base": base,
                    "Predicted State": state
                })

            df_viterbi = pd.DataFrame(viterbi_data)
            st.dataframe(df_viterbi, use_container_width=False)
            output_txt.write(df_viterbi.to_string(index=False) + '\n\n')

           # === Coding Region as Table ===
            st.subheader("=== Coding Region ===")
            output_txt.write("\n=== Coding Region ===\n")
            exon_states = ["A_start", "T_start", "G_start", "coding", "stop1", "stop2", "stop3"]
            coding_regions = []
            current = []
            region_idx = 1
            start_idx = None

            for i, state in enumerate(path_state):
                if state in exon_states:
                    if start_idx is None:
                        start_idx = i
                    current.append(obs_raw[i])
                else:
                    if current:
                        end_idx = i - 1
                        sequence = ''.join(current)
                        region_length = len(current)
                        if region_length >= min_coding_length:
                            region_info = {
                                "Start": start_idx,
                                "End": end_idx,
                                "Sequence": sequence,
                                "Length": region_length
                            }
                            coding_regions.append(region_info)
                            
                            line = f"Number[{region_idx}] Start[{start_idx}] End[{end_idx+1}] Length[{region_length}]: {sequence}"
                            region_idx +=1
                            st.write(line)
                            output_txt.write(line + "\n\n")
                            
                        
                        current = []
                        start_idx = None

            # Jika masih ada sisa coding region di akhir
            if current:
                end_idx = len(path_state) - 1
                sequence = ''.join(current)
                region_length = len(current)
    
                if region_length >= min_coding_length:
                    region_info = {
                        "Start": start_idx,
                        "End": end_idx,
                        "Sequence": sequence,
                        "Length": region_length
                    }
                    coding_regions.append(region_info)
                    
                    line = f"Number[{region_idx}] Start[{start_idx}] End[{end_idx}] Length[{region_length}]: {sequence}"
                    st.write(line)
                    output_txt.write(line + "\n\n")
                    region_idx +=1

            # df_coding = pd.DataFrame(coding_regions)
            # st.dataframe(df_coding, use_container_width=False)

            print_codon_translation8(obs_raw, path_state, Codon_DNA, output_txt)
            st.download_button(
                label="📥 Download TXT Output",
                data=output_txt.getvalue(),
                file_name="prokaryotic_analysis.txt",
                mime="text/plain"
            )

        def get_reverse_complement(seq):
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            return ''.join(complement[base] for base in reversed(seq))

        orf_frame = st.radio("🔍 Select ORF Frame", ["+1", "+2", "+3", "-1", "-2", "-3"], horizontal=True)
        if st.button("Run with ORF Frame (12 States)"):
            output_txt = StringIO()
            genome_sequence = ""
            for line in uploaded_file:
                line = line.decode("utf-8").strip()
                if not line.startswith(">"):
                    genome_sequence += line
                else:
                    desc = line   
            st.subheader("📄 File Information")
            output_txt.write("File Information\n")
            st.write(desc)
            output_txt.write(desc + '\n')
            st.write(f"Total Sequence Length: {len(genome_sequence)} bases")
            output_txt.write(f"Total Sequence Length: {len(genome_sequence)} bases\n")

            # Apply ORF frame adjustment
            if orf_frame.startswith("-"):
                genome_sequence = get_reverse_complement(genome_sequence)
            
            shift = int(orf_frame[-1]) - 1
            genome_sequence = genome_sequence[shift:]  # Apply +1, +2, +3 or reverse -1, -2, -3

            st.write("**🧬Processed Genome Sequence (Frame Applied):**")
            output_txt.write("\nProcessed Genome Sequence (Frame Applied):\n")
            st.write(genome_sequence)
            output_txt.write(genome_sequence)

            st.write("**🧪 Sequence Statistic:**")
            output_txt.write("\n\nSequence Statistic:\n")
            base_count, base_percent, df = count_each_nitrogen_base(genome_sequence)
            st.dataframe(df, use_container_width=False)
            output_txt.write(df.to_string(index=False) + '\n\n')

            # === 1. Definisi State dan Observasi ===
            states = [
                "A_start", "T_start", "G_start", "coding1", "coding2", "coding3",
                "stop1", "stop2", "stop3", "non coding1", "non coding2", "non coding3"
            ]
            n_states = len(states)
            observations = ['A', 'T', 'C', 'G']
            dna_map = {b: i for i, b in enumerate(observations)}

            # === 2. Matriks Transisi dan Emisi ===
            A = np.array([
                [0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0], #1
                [0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0], #2
                [0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0], #3
                [0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0], #4
                [0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0], #5
                [0,   0,   0, 0.8,   0,   0, 0.2,   0,   0,   0,   0,   0], #6
                [0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0], #7
                [0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0], #8
                [0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0], #9
                [0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0], #10
                [0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1], #11
                [0.2, 0,   0,   0,   0,   0,   0,   0,   0, 0.8,   0,   0] #12

            ])

            B = np.array([
                [1, 0, 0, 0], #1
                [0, 1, 0, 0], #2
                [0, 0, 0, 1], #3
                [0.25, 0.25, 0.25, 0.25], #4
                [0.25, 0.25, 0.25, 0.25], #5
                [0.25, 0.25, 0.25, 0.25], #6
                [0, 1, 0, 0], #7
                [0.67, 0, 0, 0.33], #8
                [0.67, 0, 0, 0.33], #9
                [0.25, 0.25, 0.25, 0.25], #10
                [0.25, 0.25, 0.25, 0.25], #11
                [0.25, 0.25, 0.25, 0.25], #12
            ])

            pi = np.zeros(n_states)
            pi[0] = 0.1
            pi[9] = 0.9

            # === Observasi ===
            obs_raw = genome_sequence[:dat_length]
            obs_seq = np.array([dna_map[b] for b in obs_raw])

            # === Viterbi ===
            path_idx = viterbi(obs_seq, A, B, pi)
            path_state = [states[i] for i in path_idx]

            # Output Matrix dan Result
            st.subheader("=== Matriks Transisi (A) ===")
            output_txt.write("\nMatriks Transisi (A)\n")
            df_A = pd.DataFrame(A, columns=states, index=states)
            st.dataframe(df_A, use_container_width=False)
            output_txt.write(df_A.to_string(index=False) + '\n\n')

            st.subheader("=== Matriks Emisi (B) ===")
            output_txt.write("\nMatriks Emisi (B)\n")
            df_B = pd.DataFrame(B, columns=observations, index=states)
            st.dataframe(df_B, use_container_width=False)
            output_txt.write(df_B.to_string(index=False) + '\n\n')

            st.subheader("=== Probabilitas Awal (π) ===")
            output_txt.write("\nProbabilitas Awal (π)\n")
            df_pi = pd.DataFrame({
                "State": states,
                "Initial Probability": pi
            })
            st.dataframe(df_pi, use_container_width=False)
            output_txt.write(df_pi.to_string(index=False) + '\n\n')

            # Hidden State Result
            st.subheader("=== Hidden State Result ===")
            output_txt.write("\nHidden State Result\n")
            viterbi_data = []
            for i, (base, state) in enumerate(zip(obs_raw, path_state)):
                viterbi_data.append({
                    "Position": i,
                    "Base": base,
                    "Predicted State": state
                })

            df_viterbi = pd.DataFrame(viterbi_data)
            st.dataframe(df_viterbi, use_container_width=False)
            output_txt.write(df_viterbi.to_string(index=False) + '\n\n')

            # Coding Region Table
            st.subheader("=== Coding Region ===")
            output_txt.write("\n=== Coding Region ===\n")
            exon_states = ["A_start", "T_start", "G_start", "coding1","coding2","coding3", "stop1", "stop2", "stop3"]
            coding_regions = []
            current = []
            region_idx = 1
            start_idx = None

            for i, state in enumerate(path_state):
                if state in exon_states:
                    if start_idx is None:
                        start_idx = i
                    current.append(obs_raw[i])
                else:
                    if current:
                        end_idx = i - 1
                        if len(current)> min_coding_length:
                            sequence = ''.join(current)
                            line = f"Number[{region_idx}] Start[{start_idx}] End[{end_idx}] Length[{len(current)}]: {sequence}"
                            region_idx+=1
                            st.write(line)
                            output_txt.write(line + "\n")
                        current = []
                        start_idx = None

            if current:
                end_idx = len(path_state) - 1
                if len(current)> min_coding_length:
                    sequence = ''.join(current)
                    line = f"Number[{region_idx}] Start[{start_idx}] End[{end_idx}] Length[{len(current)}]: {sequence}"
                    region_idx+=1
                    st.write(line)
                    output_txt.write(line + "\n")

            print_codon_translation12(obs_raw, path_state, Codon_DNA, output_txt)

            st.download_button(
                label="📥 Download TXT Output",
                data=output_txt.getvalue(),
                file_name="defined_sequence_analysis.txt",
                mime="text/plain"
            )
            
        # if st.button("Run with Defined Sequence"):
        #     output_txt = StringIO()
        #     desc = "Predefined Genome Sequence"
        #     genome_sequence = ['C','C','G','G','C','A','C','T','G','T','T','C','A','T','G','G','G','C','A','A',
        #            'T','G','C','A','A','G','G','T','A','C','G','G','T','G','A','G','C','A','G','G',
        #            'T','A','A','G','T','G','A','T','T','A','A','T','G','C','A','T','T','T','C','T',
        #            'C','G','C','C','A','G','T','G','G','C','T','A','G','A','C','G','A','T','G','C',
        #            'A','T','A','G','G','A','G','A','T','C','A','T','T','G','A','C','G','A','T','G',
        #            'C','A','T','A','G','A','C','C','G','G','A','A','G','C']

        #     # Join list ke string
        #     genome_sequence_str = ''.join(genome_sequence)

        #     st.subheader("📄 File Information")
        #     output_txt.write("File Information\n")
        #     st.write(desc)
        #     output_txt.write(desc + '\n')
        #     st.write(f"Total Sequence Length: {len(genome_sequence)} bases")
        #     output_txt.write(f"Total Sequence Length: {len(genome_sequence)} bases\n")

        #     st.write("**🧬Complete Genome Sequence:**")
        #     output_txt.write("\nComplete Genome Sequence:\n")
        #     st.write(genome_sequence_str)
        #     output_txt.write(genome_sequence_str)

        #     st.write("**🧪 Sequence Statistic:**")
        #     output_txt.write("\n\nSequence Statistic:\n")
        #     base_count, base_percent, df = count_each_nitrogen_base(genome_sequence)
        #     st.dataframe(df, use_container_width=False)
        #     output_txt.write(df.to_string(index=False) + '\n\n')

        #     # === 1. Definisi State dan Observasi ===
        #     states = [
        #         "A_start", "T_start", "G_start", "coding",
        #         "stop1", "stop2", "stop3", "non coding"
        #     ]
        #     n_states = len(states)
        #     observations = ['A', 'T', 'C', 'G']
        #     dna_map = {b: i for i, b in enumerate(observations)}

        #     # === 2. Matriks Transisi dan Emisi ===
        #     A = np.array([
        #         [0,   1,   0,   0,   0,   0,   0,   0],
        #         [0,   0,   1,   0,   0,   0,   0,   0],
        #         [0,   0,   0,   1,   0,   0,   0,   0],
        #         [0,   0,   0, 0.9, 0.1,   0,   0,   0],
        #         [0,   0,   0,   0,   0,   1,   0,   0],
        #         [0,   0,   0,   0,   0,   0,   1,   0],
        #         [0,   0,   0,   0,   0,   0,   0,   1],
        #         [0.2, 0,   0,   0,   0,   0,   0, 0.8]
        #     ])

        #     B = np.array([
        #         [1, 0, 0, 0],
        #         [0, 1, 0, 0],
        #         [0, 0, 0, 1],
        #         [0.25, 0.25, 0.25, 0.25],
        #         [0, 1, 0, 0],
        #         [0.67, 0, 0, 0.33],
        #         [0.67, 0, 0, 0.33],
        #         [0.25, 0.25, 0.25, 0.25],
        #     ])

        #     pi = np.zeros(n_states)
        #     pi[0] = 0.1
        #     pi[7] = 0.9

        #     # === Observasi ===
        #     obs_raw = genome_sequence[:dat_length]
        #     obs_seq = np.array([dna_map[b] for b in obs_raw])

        #     # === Viterbi ===
        #     path_idx = viterbi(obs_seq, A, B, pi)
        #     path_state = [states[i] for i in path_idx]

        #     # Output Matrix dan Result
        #     st.subheader("=== Matriks Transisi (A) ===")
        #     output_txt.write("\nMatriks Transisi (A)\n")
        #     df_A = pd.DataFrame(A, columns=states, index=states)
        #     st.dataframe(df_A, use_container_width=False)
        #     output_txt.write(df_A.to_string(index=False) + '\n\n')

        #     st.subheader("=== Matriks Emisi (B) ===")
        #     output_txt.write("\nMatriks Emisi (B)\n")
        #     df_B = pd.DataFrame(B, columns=observations, index=states)
        #     st.dataframe(df_B, use_container_width=False)
        #     output_txt.write(df_B.to_string(index=False) + '\n\n')

        #     st.subheader("=== Probabilitas Awal (π) ===")
        #     output_txt.write("\nProbabilitas Awal (π)\n")
        #     df_pi = pd.DataFrame({
        #         "State": states,
        #         "Initial Probability": pi
        #     })
        #     st.dataframe(df_pi, use_container_width=False)
        #     output_txt.write(df_pi.to_string(index=False) + '\n\n')

        #     # Hidden State Result
        #     st.subheader("=== Hidden State Result ===")
        #     output_txt.write("\nHidden State Result\n")
        #     viterbi_data = []
        #     for i, (base, state) in enumerate(zip(obs_raw, path_state)):
        #         viterbi_data.append({
        #             "Position": i,
        #             "Base": base,
        #             "Predicted State": state
        #         })

        #     df_viterbi = pd.DataFrame(viterbi_data)
        #     st.dataframe(df_viterbi, use_container_width=False)
        #     output_txt.write(df_viterbi.to_string(index=False) + '\n\n')

        #     # Coding Region Table
        #     st.subheader("=== Coding Region ===")
        #     output_txt.write("\n=== Coding Region ===\n")
        #     exon_states = ["A_start", "T_start", "G_start", "coding", "stop1", "stop2", "stop3"]
        #     coding_regions = []
        #     current = []
        #     start_idx = None

        #     for i, state in enumerate(path_state):
        #         if state in exon_states:
        #             if start_idx is None:
        #                 start_idx = i
        #             current.append(obs_raw[i])
        #         else:
        #             if current:
        #                 end_idx = i - 1
        #                 sequence = ''.join(current)
        #                 line = f"Start[{start_idx}] End[{end_idx}] Length[{len(current)}]: {sequence}"
        #                 st.write(line)
        #                 output_txt.write(line + "\n")
        #                 current = []
        #                 start_idx = None

        #     if current:
        #         end_idx = len(path_state) - 1
        #         sequence = ''.join(current)
        #         line = f"Start[{start_idx}] End[{end_idx}] Length[{len(current)}]: {sequence}"
        #         st.write(line)
        #         output_txt.write(line + "\n")

        #     print_codon_translation8(obs_raw, path_state, Codon_DNA, output_txt)

        #     st.download_button(
        #         label="📥 Download TXT Output",
        #         data=output_txt.getvalue(),
        #         file_name="defined_sequence_analysis.txt",
        #         mime="text/plain"
        #     )

        # if st.button("Run with Defined Sequence (Modif)"):
        #     output_txt = StringIO()
        #     desc = "Predefined Genome Sequence"
        #     genome_sequence = ['C','C','G','G','C','A','C','T','G','T','T','C','A','T','G','G','G','C','A','A',
        #            'T','G','C','A','A','G','G','T','A','C','G','G','T','G','A','G','C','A','G','G',
        #            'T','A','A','G','T','G','A','T','T','A','A','T','G','C','A','T','T','T','C','T',
        #            'C','G','C','C','A','G','T','G','G','C','T','A','G','A','C','G','A','T','G','C',
        #            'A','T','A','G','G','A','G','A','T','C','A','T','T','G','A','C','G','A','T','G',
        #            'C','A','T','A','G','A','C','C','G','G','A','A','G','C']

        #     # Join list ke string
        #     genome_sequence_str = ''.join(genome_sequence)

        #     st.subheader("📄 File Information")
        #     output_txt.write("File Information\n")
        #     st.write(desc)
        #     output_txt.write(desc + '\n')
        #     st.write(f"Total Sequence Length: {len(genome_sequence)} bases")
        #     output_txt.write(f"Total Sequence Length: {len(genome_sequence)} bases\n")

        #     st.write("**🧬Complete Genome Sequence:**")
        #     output_txt.write("\nComplete Genome Sequence:\n")
        #     st.write(genome_sequence_str)
        #     output_txt.write(genome_sequence_str)

        #     st.write("**🧪 Sequence Statistic:**")
        #     output_txt.write("\n\nSequence Statistic:\n")
        #     base_count, base_percent, df = count_each_nitrogen_base(genome_sequence)
        #     st.dataframe(df, use_container_width=False)
        #     output_txt.write(df.to_string(index=False) + '\n\n')

        #     # === 1. Definisi State dan Observasi ===
        #     states = [
        #         "A_start", "T_start", "G_start", "coding1", "coding2", "coding3",
        #         "stop1", "stop2", "stop3", "non coding1", "non coding2", "non coding3"
        #     ]
        #     n_states = len(states)
        #     observations = ['A', 'T', 'C', 'G']
        #     dna_map = {b: i for i, b in enumerate(observations)}

        #     # === 2. Matriks Transisi dan Emisi ===
        #     A = np.array([
        #         [0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0], #1
        #         [0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0], #2
        #         [0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0], #3
        #         [0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0], #4
        #         [0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0], #5
        #         [0,   0,   0, 0.8,   0,   0, 0.2,   0,   0,   0,   0,   0], #6
        #         [0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0], #7
        #         [0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0], #8
        #         [0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0], #9
        #         [0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0], #10
        #         [0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1], #11
        #         [0.2, 0,   0,   0,   0,   0,   0,   0,   0, 0.8,   0,   0] #12

        #     ])

        #     B = np.array([
        #         [1, 0, 0, 0], #1
        #         [0, 1, 0, 0], #2
        #         [0, 0, 0, 1], #3
        #         [0.25, 0.25, 0.25, 0.25], #4
        #         [0.25, 0.25, 0.25, 0.25], #5
        #         [0.25, 0.25, 0.25, 0.25], #6
        #         [0, 1, 0, 0], #7
        #         [0.67, 0, 0, 0.33], #8
        #         [0.67, 0, 0, 0.33], #9
        #         [0.25, 0.25, 0.25, 0.25], #10
        #         [0.25, 0.25, 0.25, 0.25], #11
        #         [0.25, 0.25, 0.25, 0.25], #12
        #     ])

        #     pi = np.zeros(n_states)
        #     pi[0] = 0.1
        #     pi[9] = 0.9

        #     # === Observasi ===
        #     obs_raw = genome_sequence[:dat_length]
        #     obs_seq = np.array([dna_map[b] for b in obs_raw])

        #     # === Viterbi ===
        #     path_idx = viterbi(obs_seq, A, B, pi)
        #     path_state = [states[i] for i in path_idx]

        #     # Output Matrix dan Result
        #     st.subheader("=== Matriks Transisi (A) ===")
        #     output_txt.write("\nMatriks Transisi (A)\n")
        #     df_A = pd.DataFrame(A, columns=states, index=states)
        #     st.dataframe(df_A, use_container_width=False)
        #     output_txt.write(df_A.to_string(index=False) + '\n\n')

        #     st.subheader("=== Matriks Emisi (B) ===")
        #     output_txt.write("\nMatriks Emisi (B)\n")
        #     df_B = pd.DataFrame(B, columns=observations, index=states)
        #     st.dataframe(df_B, use_container_width=False)
        #     output_txt.write(df_B.to_string(index=False) + '\n\n')

        #     st.subheader("=== Probabilitas Awal (π) ===")
        #     output_txt.write("\nProbabilitas Awal (π)\n")
        #     df_pi = pd.DataFrame({
        #         "State": states,
        #         "Initial Probability": pi
        #     })
        #     st.dataframe(df_pi, use_container_width=False)
        #     output_txt.write(df_pi.to_string(index=False) + '\n\n')

        #     # Hidden State Result
        #     st.subheader("=== Hidden State Result ===")
        #     output_txt.write("\nHidden State Result\n")
        #     viterbi_data = []
        #     for i, (base, state) in enumerate(zip(obs_raw, path_state)):
        #         viterbi_data.append({
        #             "Position": i,
        #             "Base": base,
        #             "Predicted State": state
        #         })

        #     df_viterbi = pd.DataFrame(viterbi_data)
        #     st.dataframe(df_viterbi, use_container_width=False)
        #     output_txt.write(df_viterbi.to_string(index=False) + '\n\n')

        #     # Coding Region Table
        #     st.subheader("=== Coding Region ===")
        #     output_txt.write("\n=== Coding Region ===\n")
        #     exon_states = ["A_start", "T_start", "G_start", "coding1","coding2","coding3", "stop1", "stop2", "stop3"]
        #     coding_regions = []
        #     current = []
        #     start_idx = None

        #     for i, state in enumerate(path_state):
        #         if state in exon_states:
        #             if start_idx is None:
        #                 start_idx = i
        #             current.append(obs_raw[i])
        #         else:
        #             if current:
        #                 end_idx = i - 1
        #                 sequence = ''.join(current)
        #                 line = f"Start[{start_idx}] End[{end_idx}] Length[{len(current)}]: {sequence}"
        #                 st.write(line)
        #                 output_txt.write(line + "\n")
        #                 current = []
        #                 start_idx = None

        #     if current:
        #         end_idx = len(path_state) - 1
        #         sequence = ''.join(current)
        #         line = f"Start[{start_idx}] End[{end_idx}] Length[{len(current)}]: {sequence}"
        #         st.write(line)
        #         output_txt.write(line + "\n")

        #     print_codon_translation12(obs_raw, path_state, Codon_DNA, output_txt)

        #     st.download_button(
        #         label="📥 Download TXT Output",
        #         data=output_txt.getvalue(),
        #         file_name="defined_sequence_analysis.txt",
        #         mime="text/plain"
        #     )

        if st.button("Run Eukaryotic"):
            output_txt = StringIO()
            genome_sequence = ['C','C','G','G','C','A','C','T','G','T','T','C','A','T','G','G','G','C','A','A',
                    'T','G','C','A','A','G','G','T','A','C','G','G','T','G','A','G','C','A','G','G',
                    'T','A','A','G','T','G','A','T','T','A','A','T','G','C','A','T','T','T','C','T',
                    'C','G','C','C','A','G','T','G','G','C','T','A','G','A','C','G','A','T','G','C',
                    'A','T','A','G','G','A','G','A','T','C','A','T','T','G','A','C','G','A','T','G',
                    'C','A','T','A','G','A','C','C','G','G','A','A','G','C']
            
            # Terapkan ORF frame
            if orf_frame.startswith("-"):
                genome_sequence = get_reverse_complement(genome_sequence)

            shift = int(orf_frame[-1]) - 1
            genome_sequence = genome_sequence[shift:]
    
            st.subheader("📄 File Information")
            output_txt.write("File Information\n")
            st.write(f"Total Sequence Length: {len(genome_sequence)} bases")
            output_txt.write(f"Total Sequence Length: {len(genome_sequence)} bases\n")

            st.write("**🧬Processed Genome Sequence (Frame Applied):**")
            output_txt.write("\nProcessed Genome Sequence (Frame Applied):\n")
            st.write("".join(genome_sequence))
            output_txt.write("".join(genome_sequence))

            st.write("**🧪 Sequence Statistic:**")
            output_txt.write("\n\nSequence Statistic:\n")
            base_count, base_percent, df = count_each_nitrogen_base(genome_sequence)
            st.dataframe(df, use_container_width=False)
            output_txt.write(df.to_string(index=False) + '\n\n')

            # === 1. Definisi State dan Observasi ===
            states = [
                "A_start", "T_start", "G_start",
                "Coding1", "Coding2", "Coding3",
                "Stop1", "Stop2", "Stop3",
                "NC1", "NC2", "NC3",
                "G_NC", "T_NC", "A_NC",
                "C_NC2", "A_NC2", "G_NC2"
            ]
            n_states = len(states)
            observations = ['A', 'T', 'C', 'G']
            dna_map = {b: i for i, b in enumerate(observations)}

            # === 2. Matriks Transisi ===
            A = np.array([
                [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #A_start
                [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #T_start
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #G_start
                [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #code1
                [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #code2
                [0, 0, 0, 0.25, 0, 0, 0.4, 0, 0, 0, 0, 0, 0.35, 0, 0, 0, 0, 0], #code3
                [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #stop1
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], #stop2
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], #stop3
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], #nc1
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], #nc2
                [0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0, 0.1, 0, 0], #nc3
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], #g nc
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], #t nc
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], #a nc
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], #c nc2
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], #a nc2
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #g nc2
            ])

            # === 3. Matriks Emisi ===
            B = np.array([
                [1, 0, 0, 0],      # A_start
                [0, 1, 0, 0],      # T_start
                [0, 0, 0, 1],      # G_start
                [0.25, 0.25, 0.25, 0.25],  # Coding1
                [0.25, 0.25, 0.25, 0.25],  # Coding2
                [0.25, 0.25, 0.25, 0.25],  # Coding3
                [0, 1, 0, 0],      # Stop1
                [0.67, 0, 0, 0.33],# Stop2
                [0.67, 0, 0, 0.33],# Stop3
                [0.25, 0.25, 0.25, 0.25],  # NC1
                [0.25, 0.25, 0.25, 0.25],  # NC2
                [0.25, 0.25, 0.25, 0.25],  # NC3
                [0, 0, 0, 1],      # G_NC
                [0, 1, 0, 0],      # T_NC
                [1, 0, 0, 0],      # A_NC
                [0, 0, 1, 0],      # C_NC2
                [1, 0, 0, 0],      # A_NC2
                [0, 0, 0, 1],      # G_NC2
            ])


            # === 4. Probabilitas Awal ===
            pi = np.zeros(n_states)
            pi[0] = 0.1  # A_start
            pi[9] = 0.9  # NC1

            # === 5. Observasi dari user ===
 
            obs_raw = genome_sequence
            obs_seq = np.array([dna_map[b] for b in obs_raw])

            # === 6. Run Viterbi ===
            path_idx = viterbi(obs_seq, A, B, pi)
            path_state = [states[i] for i in path_idx]


            # === 7. Print Sample Output ===
            # === Print Matriks Transisi A ===
            st.subheader("=== Matriks Transisi (A) ===")
            output_txt.write("\nMatriks Transisi (A)\n")
            df_A = pd.DataFrame(A, columns=states, index=states)
            st.dataframe(df_A, use_container_width=False)
            output_txt.write(df_A.to_string(index=False) + '\n\n')

            # === Print Matriks Emisi B ===
            st.subheader("=== Matriks Emisi (B) ===")
            output_txt.write("\nMatriks Emisi (B)\n")
            df_B = pd.DataFrame(B, columns=["A", "T", "C", "G"], index=states)
            st.dataframe(df_B, use_container_width=False)
            output_txt.write(df_B.to_string(index=False) + '\n\n')


            # === Print Probabilitas Awal π ===
            st.subheader("=== Probabilitas Awal (π) ===")
            output_txt.write("\nProbabilitas Awal (π)\n")
            df_pi = pd.DataFrame({
                "State": states,
                "Initial Probability": pi
            })
            st.dataframe(df_pi, use_container_width=False)
            output_txt.write(df_pi.to_string(index=False) + '\n\n')


            # === Print Semua State dari Viterbi ===
            st.subheader("=== Hidden State Result ===")
            output_txt.write("\nHidden State Result\n")
            current = []
            viterbi_data = []
            for i, (base, state) in enumerate(zip(obs_raw, path_state)):
                viterbi_data.append({
                    "Position": i,
                    "Base": base,
                    "Predicted State": state
                })

            df_viterbi = pd.DataFrame(viterbi_data)
            st.dataframe(df_viterbi, use_container_width=False)
            output_txt.write(df_viterbi.to_string(index=False) + '\n\n')

            # === Print Exon dan Index-nya ===
            st.subheader("\n=== Segment Exon ===")
            output_txt.write("\nSegment Exon\n")
            exon_states = ["A_start", "T_start", "G_start", "Coding1", "Coding2", "Coding3", "Stop1", "Stop2", "Stop3"]
            exon_segments = []
            # exon_ranges = []
            current = []
            start_idx = None
            for i, state in enumerate(path_state):
                if state in exon_states:
                    if start_idx is None:
                        start_idx = i
                    current.append(obs_raw[i])
                else:
                    if current:
                        exon_segments.append({
                            "Start Index": start_idx-1,
                            "End Index": i - 1,
                            "Sequence": ''.join(current),
                            "Length": len(current)
                        })
                        # exon_ranges.append((start_idx, i - 1))
                        # Print per region
                        line = f"Start[{start_idx}] End[{i-1}] Length[{len(current)}]: {''.join(current)}"
                        st.write(line)
                        output_txt.write(line + "\n")

                        current = []
                        start_idx = None

            # Jika ada sisa exon di akhir
            if current:
                exon_segments.append({
                    "Start Index": start_idx-1,
                    "End Index": len(path_state) - 1,
                    "Sequence": ''.join(current),
                    "Length": len(current)
                })
                # exon_ranges.append((start_idx, len(path_state) - 1))
                # Print terakhir
                line = f"Start[{start_idx}] End[{len(path_state)-1}] Length[{len(current)}]: {''.join(current)}"
                st.write(line)
                output_txt.write(line + "\n")
            # if exon_segments:
            #     df = pd.DataFrame(exon_segments)
            #     st.dataframe(df)
            # else:
            #     st.warning("Tidak ditemukan segmen exon.")


            print_codon_translation(obs_raw, path_state, Codon_DNA, output_txt)
            st.download_button(
                label="📥 Download TXT Output",
                data=output_txt.getvalue(),
                file_name="eukaryotic_analysis.txt",
                mime="text/plain"
            )
            





    
if __name__ == "__main__":
    main()

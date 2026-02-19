# Import packages
import random
import numpy as np
import bamnostic as bs
import seqlogo

# Import modules
from data_readers import get_gff, get_fasta
from seq_ops import reverse_complement, get_seq
from motif_ops import build_pfm, build_pwm, score_kmer, pfm_ic


def initialize_random_motifs(seqs, k, rng):
    """
    Initialize random motif instances from each sequence.

    Args:
        seqs (list of str): DNA sequences
        k (int): motif length
        rng: numpy random number generator

    Returns:
        list of str: randomly selected k-mers from each sequence
    """
    # Initialize empty list to store motif instances
    motifs = []

    # Iterate through each sequence
    for seq in seqs:
        # Skip sequences shorter than motif length
        if len(seq) < k:
            continue

        # Calculate maximum valid starting position
        max_start = len(seq) - k

        # Randomly select a starting position
        start = rng.integers(0, max_start + 1)

        # Extract k-mer and add to motifs list
        motifs.append(seq[start:start + k])

    return motifs


def score_all_kmers(seq, k, pwm, rng):
    """
    Score all possible k-mers in a sequence on both strands.

    Args:
        seq (str): DNA sequence
        k (int): motif length
        pwm (np.ndarray): Position Weight Matrix
        rng: numpy random number generator

    Returns:
        tuple: (scores, strands) - lists of chosen scores and strand indicators
    """
    # Initialize lists to store scores and strand choices
    scores = []
    strands = []

    # Get sequence length
    L = len(seq)

    # Return empty lists if sequence is too short
    if L < k:
        return scores, strands

    # Iterate through all possible k-mer positions
    for start in range(0, L - k + 1):
        # Extract k-mer at current position
        kmer = seq[start:start + k]

        # Get reverse complement of k-mer
        rc = reverse_complement(kmer)

        # Score forward strand k-mer
        try:
            score_fwd = score_kmer(kmer, pwm)
        except ValueError:
            score_fwd = -np.inf

        # Score reverse strand k-mer
        try:
            score_rev = score_kmer(rc, pwm)
        except ValueError:
            score_rev = -np.inf

        # Convert scores to probabilities using softmax
        raw = np.array([score_fwd, score_rev], dtype=float)
        raw -= np.max(raw)  # Subtract max for numerical stability
        probs = np.exp(raw)
        probs /= probs.sum()  # Normalize to sum to 1

        # Randomly sample strand based on probabilities (0=forward, 1=reverse)
        strand_choice = rng.choice([0, 1], p=probs)

        # Store the score for the chosen strand
        chosen_score = score_fwd if strand_choice == 0 else score_rev
        scores.append(chosen_score)
        strands.append(strand_choice)

    return scores, strands


def scores_to_probabilities(scores):
    """
    Convert scores to probability distribution using softmax.

    Args:
        scores (list or array): numerical scores

    Returns:
        np.ndarray: probability distribution
    """
    # Convert to numpy array
    scores = np.array(scores, dtype=float)

    # Return empty array if no scores
    if scores.size == 0:
        return scores

    # Apply softmax transformation for numerical stability
    scores -= np.max(scores)
    exp_scores = np.exp(scores)

    # Calculate sum for normalization
    total = np.sum(exp_scores)

    # Handle case where all scores are very negative
    if total == 0:
        return np.ones_like(exp_scores) / len(exp_scores)

    # Normalize to create probability distribution
    return exp_scores / total


def GibbsMotifFinder(seqs, k, seed=None, max_iters=1000):
    """
    Gibbs sampler for motif convergence.

    Args:
        seqs (list of str): DNA sequences
        k (int): motif length
        seed (int): random seed
        max_iters (int): number of iterations

    Returns:
        np.ndarray: final Position Frequency Matrix (4 x k)
    """
    # Set random seed for reproducibility
    random.seed(seed)
    rng = np.random.default_rng(seed)

    # Initialize random motif instances from each sequence
    motifs = initialize_random_motifs(seqs, k, rng)

    # Validate that motifs were successfully initialized
    if len(motifs) == 0:
        raise ValueError("No motifs initialized. Check sequence lengths and k.")

    # Run Gibbs sampling iterations
    for it in range(max_iters):
        # Get number of motif instances
        N = len(motifs)

        # Randomly select one sequence to update
        i = rng.integers(0, N)

        # Build PFM from all motifs except the selected one
        motifs_except_i = [motifs[j] for j in range(N) if j != i]
        pfm = build_pfm(motifs_except_i, k)

        # Build PWM from the PFM
        pwm = build_pwm(pfm)

        # Score all k-mers in the selected sequence
        seq_i = seqs[i]
        scores, strands = score_all_kmers(seq_i, k, pwm, rng)

        # Skip if no valid k-mers were found
        if not scores:
            continue

        # Convert scores to probability distribution
        probs = scores_to_probabilities(scores)

        # Create array of position indices
        positions = np.arange(len(scores))

        # Sample a position based on probabilities
        m = rng.choice(positions, p=probs)

        # Get the strand choice for the sampled position
        chosen_strand = strands[m]

        # Update motif instance based on chosen strand
        if chosen_strand == 0:
            # Use forward strand
            motifs[i] = seq_i[m:m + k]
        else:
            # Use reverse strand
            motifs[i] = reverse_complement(seq_i[m:m + k])

    # Build final PFM from converged motifs
    final_pfm = build_pfm(motifs, k)

    return final_pfm


def load_sequences_from_bam(bam_path, max_reads=10000):
    """
    Load up to max_reads sequences from BAM file.

    Args:
        bam_path (str): path to BAM file
        max_reads (int): maximum number of reads to load

    Returns:
        list of str: DNA sequences
    """
    # Initialize list to store sequences
    seqs = []

    # Open BAM file and extract sequences
    with bs.AlignmentFile(bam_path) as bam:
        for idx, read in enumerate(bam):
            # Stop after reaching maximum reads
            if idx >= max_reads:
                break

            # Skip reads without sequence data
            if read.seq is None:
                continue

            # Add sequence to list
            seqs.append(read.seq)

    return seqs


# Main script execution
if __name__ == "__main__":
    # Set input parameters
    bam_path = "SRR9090854.subsampled_5pct.bam"  # Input BAM file path
    k = 10  # Motif length to search for
    seed = 42  # Random seed for reproducibility

    # Alternative BAM-based sequence loading (commented out)
    # print("Loading sequences from BAM file...")
    # seqs = load_sequences_from_bam(bam_path, max_reads=500)
    # print(f"Loaded {len(seqs)} sequences.")

    # Set paths for genomic data files
    seq_file = "data/GCF_000009045.1_ASM904v1_genomic.fna"
    gff_file = "data/GCF_000009045.1_ASM904v1_genomic.gff"

    # Initialize sequence list
    seqs = []

    # Load sequences from FASTA file
    for name, seq in get_fasta(seq_file):
        # Parse GFF annotations
        for gff_entry in get_gff(gff_file):
            # Filter for coding sequences only
            if gff_entry.type == "CDS":
                # Extract promoter region (50bp upstream)
                promoter_seq = get_seq(seq, gff_entry.start, gff_entry.end, gff_entry.strand, 50)

                # Filter for sequences containing ribosome binding site motif
                if "AGGAGG" in promoter_seq:
                    seqs.append(promoter_seq)

    # Run Gibbs sampling motif finder
    print("Running GibbsMotifFinder...")
    pfm = GibbsMotifFinder(seqs, k, seed=seed, max_iters=10000)
    print("Gibbs sampling complete.")

    # Display Position Frequency Matrix
    print("\nFinal PFM (4 x k):")
    print(pfm)

    # Calculate Position Weight Matrix by normalizing columns
    pwm = pfm / pfm.sum(axis=0)  # Normalize each column to sum to 1
    print("\nFinal PWM (A,C,G,T rows):")
    print(pwm)

    # Determine consensus sequence from PWM
    bases = np.array(["A", "C", "G", "T"])
    consensus = "".join(bases[np.argmax(pwm, axis=0)])  # Select highest probability base per position
    print("\nMost probable consensus sequence:", consensus)

    # Calculate information content of the motif
    ic = pfm_ic(pfm)
    print(f"\nFinal PFM information content: {ic:.3f}")

    # Write results to output text file
    with open("p53_results.txt", "w") as f:
        # Write Position Frequency Matrix
        f.write("Final PFM (4 x k):\n")
        f.write(str(pfm) + "\n\n")

        # Write Position Weight Matrix
        f.write("Final PWM (A,C,G,T rows):\n")
        f.write(str(pwm) + "\n\n")

        # Write consensus sequence
        f.write("Consensus sequence:\n")
        f.write(consensus + "\n\n")

        # Write information content
        f.write(f"Final PFM information content: {ic:.3f}\n")

    print("Saved results to p53_results.txt")

    # Generate sequence logo visualization
    print("Generating sequence logo...")
    logo_pm = seqlogo.CompletePm(pfm=pfm.T)  # Transpose PFM for seqlogo format
    seqlogo.seqlogo(logo_pm, format='png', filename='p53_logo.png')
    print("Saved sequence logo to p53_logo.png")

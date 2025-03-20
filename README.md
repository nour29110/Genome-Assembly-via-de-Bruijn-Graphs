# Implmenting an Genome Assembler: Phi-X174 genome using Overlap Graph, Kmer Composition and De-Bruijn Graph.

Phi X 174 (ΦX174) is a single-stranded DNA virus that holds the distinction of being the first DNA-based genome to be sequenced, a milestone achieved by Fred Sanger’s team in 1977. Prior to this, in 1962, Walter Fiers and Robert Sinsheimer had already confirmed that ΦX174’s DNA is a covalently closed circular molecule. Later, Nobel laureate Arthur Kornberg used ΦX174 as a model system to demonstrate that DNA synthesized in vitro by purified enzymes could recreate all the attributes of a natural virus—paving the way for synthetic biology. In 2003, Craig Venter’s group successfully assembled the entire ΦX174 genome in vitro from synthetic oligonucleotides, even achieving the complete in vitro assembly of the virus particle. More recent research has shown that its highly overlapping genome can be fully decompressed while retaining full functionality.

<img width="259" alt="default" src="https://github.com/nour29110/Genome-Assembly-via-de-Bruijn-Graphs/blob/main/week3_Assembly%20Faces%20Real%20Sequencing%20Data/Phi-X174.png">

---

## Problem 1: Finding a Circulation in a Network

### Overview
Determine if a valid circulation exists in a directed network where each edge has a lower bound and capacity. The flow on each edge must:
- Satisfy the capacity constraints: _l ≤ f ≤ c_.
- Preserve flow conservation at every vertex (incoming flow equals outgoing flow).

### Key Tasks
- Model the network with vertices and directed edges.
- Assign flows respecting the given lower bounds and capacities.
- Verify flow conservation at each vertex.
- Output "YES" along with the flow values for each edge if a circulation exists; otherwise, output "NO".

---

## Problem 2: Selecting the Optimal k-mer Size

### Overview
Select an optimal k-mer length for fragmenting error-free reads such that the de Bruijn graph contains a single Eulerian cycle. The challenge lies in balancing:
- A k-mer that is small enough to ensure connectivity.
- A k-mer that is large enough to avoid an overly complex graph.

### Key Tasks
- Experiment with various k-mer sizes.
- Fragment the reads into segments of length _k_.
- Construct the corresponding de Bruijn graph.
- Identify the optimal _k_ that produces a single Eulerian cycle.
- Output the optimal _k_ value.

---

## Problem 3: Bubble Detection

### Overview
Detect bubbles in a de Bruijn graph constructed from error-prone reads. Bubbles are formed by two disjoint, non-overlapping paths (from the same start to end nodes) caused by sequencing errors.

### Key Tasks
- Break error-prone reads into k-mers.
- Construct the de Bruijn graph.
- Detect and count (v, w)-bubbles using a given bubble length threshold.
- Output the number of bubbles detected.

---

## Problem 4: Tip Removal

### Overview
Remove spurious paths (tips) from the de Bruijn graph that are a result of sequencing errors. Tips are typically branches that start at nodes with no incoming edges or end at nodes with no outgoing edges.

### Key Tasks
- Construct the de Bruijn graph from 15-mers derived from the reads.
- Identify tip structures in the graph.
- Iteratively remove the tips and update the graph accordingly.
- Output the total number of edges removed during tip removal.

---

## Problem 5: Assembling the phi X174 Genome from Error-Prone Reads

### Overview
Assemble the circular genome of the phi X174 bacteriophage using error-prone reads. This involves:
- Breaking reads into 15-mers.
- Constructing a de Bruijn graph.
- Applying error correction techniques (tip removal and bubble resolution).
- Reconstructing the circular genome.

### Key Tasks
- Fragment the error-prone reads into 15-mers.
- Build and refine the de Bruijn graph using error correction.
- Reconstruct the circular genome from the corrected graph.
- Output the assembled genome in FASTA format on a single line.

---

## Problem 6: Final Project – Implementing an Assembler

### Overview
Integrate all previous components to create a complete genome assembler that accepts reads or read-pairs as input, applies error correction, and outputs assembled contigs in FASTA format.

[genome1](https://d3c33hcgiwev3.cloudfront.net/_c1223813227b2ecec3e60224e6f070e4_genome1.txt?Expires=1742601600&Signature=RSp~QY8booeetJuQMRpnwB28Pz637t-hZM6b-zkUdWW8E-So7XkE2HohOtfLi8~x4c7qkH0JUyV-9BS2KXBSN~8y-zOQGLeXMfKSIUMdpaYOxJ~WuuVRgvxrcnQwh41z53G0ymH9XyHXav3Y3a0XQnouRM3b9dFlXQjTj-MTzyc_&Key-Pair-Id=APKAJLTNE6QMUY6HBC5A)
[dataset1](https://d3c33hcgiwev3.cloudfront.net/_c1223813227b2ecec3e60224e6f070e4_dataset1.txt?Expires=1742601600&Signature=WTKAd-dbcO6rsNRn3hPtI5UKE-ePVg1KL5NGprVvhTe-Hru-cCZlGdaX1OXoyntYicuJCAs7s4W4rVRfx5ctiLYm7bpkiioLavctgN25mq3iA0XjKZ~XYrTov7m418r75y7McIsnuA5JKdtj26IbxlW1SKMrZ7qoQBArikLBAc4_&Key-Pair-Id=APKAJLTNE6QMUY6HBC5A)

### Key Tasks
- Parse input containing either single reads or read-pairs.
- Apply various error-correction techniques (e.g., tip removal, bubble detection).
- Construct the de Bruijn graph from processed reads.
- Extract contigs from the graph.
- Output the assembled contigs in FASTA format.
- Ensure the assembler works with different datasets (error-free, error-prone, and real sequencing data).

---

Each component builds on the previous ones to form a comprehensive genome assembly pipeline. 

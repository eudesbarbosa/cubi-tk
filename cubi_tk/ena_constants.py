"""Values for controlled vocabularies from ENA.

Taken from

- https://ena-docs.readthedocs.io/en/latest/submit/reads/webin-cli.html
"""

# Constants for platform definitions.

LS454 = "LS454"
ILLUMINA = "ILLUMINA"
PACBIO_SMRT = "PACBIO_SMRT"
IONTORRENT = "ION_TORRENT"
CAPILLARY = "CAPILLARY"
ONT = "OXFORD_NANOPORE"
BGISEQ = "BGISEQ"
DNBSEQ = "DNBSEQ"

#: Allowed platforms in ENA.
PLATFORMS = (LS454, ILLUMINA, PACBIO_SMRT, IONTORRENT, CAPILLARY, ONT, BGISEQ, DNBSEQ)

# Constants for platforms.

LS454_454GS = "454 GS"
LS454_454GS_240 = "454 GS 20"
LS454_454GS_FLX = "454 GS FLX"
LS454_454GS_FLX_PLUS = "454 GS FLX+"
LS454_454GS_FLX_TITANIUM = "454 GS FLX Titanium"
LS454_454GS_JUNIOR = "454 GS Junior"
ILLUMINA_HISEQX_FIVE = "HiSeq X Five"
ILLUMINA_HISEQX_TEN = "HiSeq X Ten"
ILLUMINA_GA = "Illumina Genome Analyzer"
ILLUMINA_GA2 = "Illumina Genome Analyzer II"
ILLUMINA_GA2x = "Illumina Genome Analyzer IIx"
ILLUMINA_HISCANQ = "Illumina HiScanSQ"
ILLUMINA_HISEQ_1000 = "Illumina HiSeq 1000"
ILLUMINA_HISEQ_1500 = "Illumina HiSeq 1500"
ILLUMINA_HISEQ_2000 = "Illumina HiSeq 2000"
ILLUMINA_HISEQ_2500 = "Illumina HiSeq 2500"
ILLUMINA_HISEQ_3000 = "Illumina HiSeq 3000"
ILLUMINA_HISEQ_4000 = "Illumina HiSeq 4000"
ILLUMINA_ISEQ_100 = "Illumina iSeq 100"
ILLUMINA_MISEQ = "Illumina MiSeq"
ILLUMINA_MINISEQ = "Illumina MiniSeq"
ILLUMINA_NOVASEQ_6000 = "Illumina NovaSeq 6000"
ILLUMINA_NETSEQ_500 = "NextSeq 500"
ILLUMINA_NETSEQ_550 = "NextSeq 550"
PACBIO_RS = "PacBio RS"
PACBIO_RS2 = "PacBio RS II"
PACBIO_SEQEL = "Sequel"
IONTORRENT_PGM = "Ion Torrent PGM"
IONTORRENT_PROTON = "Ion Torrent Proton"
IONTORRENT_S5 = "Ion Torrent S5"
IONTORRENT_S5XL = "Ion Torrent S5 XL"
ABI_AB3730XL = "AB 3730xL Genetic Analyzer"
ABI_AB3730 = "AB 3730 Genetic Analyzer"
ABI_AB3500XL = "AB 3500xL Genetic Analyzer"
ABI_AB3500 = "AB 3500 Genetic Analyzer"
ABI_AB3130XL = "AB 3130xL Genetic Analyzer"
ABI_AB3130 = "AB 3130 Genetic Analyzer"
ABI_AB310 = "AB 310 Genetic Analyzer"
ONT_MINION = "MinION"
ONT_GRIDION = "GridION"
ONT_PROMETHION = "PromethION"
BGI_BGISEQ500 = "BGISEQ-500"
DNB_T7 = "DNBSEQ-T7"
DNB_G400 = "DNBSEQ-G400"
DNB_G50 = "DNBSEQ-G50"
DNB_G400_FAST = "DNBSEQ-G400 FAST"
UNSPECIFIED = "unspecified"

#: Allowed values for instruments in ENA records.
INSTRUMENTS = (
    LS454_454GS,
    LS454_454GS_240,
    LS454_454GS_FLX,
    LS454_454GS_FLX_PLUS,
    LS454_454GS_FLX_TITANIUM,
    LS454_454GS_JUNIOR,
    ILLUMINA_HISEQX_FIVE,
    ILLUMINA_HISEQX_TEN,
    ILLUMINA_GA,
    ILLUMINA_GA2,
    ILLUMINA_GA2x,
    ILLUMINA_HISCANQ,
    ILLUMINA_HISEQ_1000,
    ILLUMINA_HISEQ_1500,
    ILLUMINA_HISEQ_2000,
    ILLUMINA_HISEQ_2500,
    ILLUMINA_HISEQ_3000,
    ILLUMINA_HISEQ_4000,
    ILLUMINA_ISEQ_100,
    ILLUMINA_MISEQ,
    ILLUMINA_MINISEQ,
    ILLUMINA_NOVASEQ_6000,
    ILLUMINA_NETSEQ_500,
    ILLUMINA_NETSEQ_550,
    PACBIO_RS,
    PACBIO_RS2,
    PACBIO_SEQEL,
    IONTORRENT_PGM,
    IONTORRENT_PROTON,
    IONTORRENT_S5,
    IONTORRENT_S5XL,
    ABI_AB3730XL,
    ABI_AB3730,
    ABI_AB3500XL,
    ABI_AB3500,
    ABI_AB3130XL,
    ABI_AB3130,
    ABI_AB310,
    ONT_MINION,
    ONT_GRIDION,
    ONT_PROMETHION,
    BGI_BGISEQ500,
    DNB_T7,
    DNB_G400,
    DNB_G50,
    DNB_G400_FAST,
    UNSPECIFIED,
)

# Constants for library selection.

LIBSEL_RANDOM = "RANDOM"
LIBSEL_PCR = "PCR"
LIBSEL_RANDOM_PCR = "RANDOM PCR"
LIBSEL_RT_PCR = "RT-PCR"
LIBSEL_HMPR = "HMPR"
LIBSEL_MF = "MF"
LIBSEL_REPEAT_FRACTIONATION = "repeat fractionation"
LIBSEL_SIZE_FRACTIONATION = "size fractionation"
LIBSEL_MSLL = "MSLL"
LIBSEL_CDNA = "cDNA"
LIBSEL_CDNA_RANDOM_PRIMING = "cDNA_randomPriming"
LIBSEL_CDNA_OLIGO_DT = "cDNA_oligo_dT"
LIBSEL_POLYA = "PolyA"
LIBSEL_OLIGO_DT = "Oligo-dT"
LIBSEL_INVERSE_RNA = "Inverse rRNA"
LIBSEL_INVERSE_RNA_SELECTION = "Inverse rRNA selection"
LIBSEL_CHIP = "ChIP"
LIBSEL_CHIP_SEQ = "ChIP-Seq"
LIBSEL_MNASE = "MNase"
LIBSEL_DNASE = "DNase"
LIBSEL_HYBRID_SELECTION = "Hybrid Selection"
LIBSEL_REDUCED_REPRESENTATION = "Reduced Representation"
LIBSEL_RESTRICTION_DIGEST = "Restriction Digest"
LIBSEL_5HETYHLCYTITINE_ANTIBODY = "5-methylcytidine antibody"
LIBSEL_MBD2_PROTEIN_METHYL_CPG_BINDING_DOMAIN = "MBD2 protein methyl-CpG binding domain"
LIBSEL_CAGE = "CAGE"
LIBSEL_RACE = "RACE"
LIBSEL_MDA = "MDA"
LIBSEL_PADLOCK_PROBES_CPATURE_METHOD = "padlock probes capture method"
LIBSEL_OTHER = "other"
LIBSEL_UNSPECIFIED = "unspecified"

#: Allowed library selection strategies for ENA records.
LIBRARY_SELECTIONS = (
    LIBSEL_RANDOM,
    LIBSEL_PCR,
    LIBSEL_RANDOM_PCR,
    LIBSEL_RT_PCR,
    LIBSEL_HMPR,
    LIBSEL_MF,
    LIBSEL_REPEAT_FRACTIONATION,
    LIBSEL_SIZE_FRACTIONATION,
    LIBSEL_MSLL,
    LIBSEL_CDNA,
    LIBSEL_CDNA_RANDOM_PRIMING,
    LIBSEL_CDNA_OLIGO_DT,
    LIBSEL_POLYA,
    LIBSEL_OLIGO_DT,
    LIBSEL_INVERSE_RNA,
    LIBSEL_INVERSE_RNA_SELECTION,
    LIBSEL_CHIP,
    LIBSEL_CHIP_SEQ,
    LIBSEL_MNASE,
    LIBSEL_DNASE,
    LIBSEL_HYBRID_SELECTION,
    LIBSEL_REDUCED_REPRESENTATION,
    LIBSEL_RESTRICTION_DIGEST,
    LIBSEL_5HETYHLCYTITINE_ANTIBODY,
    LIBSEL_MBD2_PROTEIN_METHYL_CPG_BINDING_DOMAIN,
    LIBSEL_CAGE,
    LIBSEL_RACE,
    LIBSEL_MDA,
    LIBSEL_PADLOCK_PROBES_CPATURE_METHOD,
    LIBSEL_OTHER,
    LIBSEL_UNSPECIFIED,
)

# Constants for library sources.

LIBSRC_GENOMIC = "GENOMIC"
LIBSRC_GENOMIC_SC = "GENOMIC SINGLE CELL"
LIBSRC_TRANSCRIPTOMIC = "TRANSCRIPTOMIC"
LIBSRC_TRANSCRIPTOMIC_SC = "TRANSCRIPTOMIC SINGLE CELL"
LIBSRC_METAGENOMIC = "METAGENOMIC"
LIBSRC_METATRANSCRIPTOMIC = "METATRANSCRIPTOMIC"
LIBSRC_SYNTHETIC = "SYNTHETIC"
LIBSRC_VIRAL_RNA = "VIRAL RNA"
LIBSRC_OTHER = "OTHER"

#: Allowed library sources for ENA records.
LIBRARY_SOURCES = (
    LIBSRC_GENOMIC,
    LIBSRC_GENOMIC_SC,
    LIBSRC_TRANSCRIPTOMIC,
    LIBSRC_TRANSCRIPTOMIC_SC,
    LIBSRC_METAGENOMIC,
    LIBSRC_METATRANSCRIPTOMIC,
    LIBSRC_SYNTHETIC,
    LIBSRC_VIRAL_RNA,
    LIBSRC_OTHER,
)

# Constants for library strategies.

LIBSTR_WGS = "WGS"
LIBSTR_WGA = "WGA"
LIBSTR_WXS = "WXS"
LIBSTR_RNA_SEQ = "RNA-Seq"
LIBSTR_SSRNA_SEQ = "ssRNA-seq"
LIBSTR_MIRNA_SEQ = "miRNA-Seq"
LIBSTR_NCRNA_SEQ = "ncRNA-Seq"
LIBSTR_FL_CDNA = "FL-cDNA"
LIBSTR_EST = "EST"
LIBSTR_HIC = "Hi-C"
LIBSTR_ATAC_SEQ = "ATAC-seq"
LIBSTR_WCS = "WCS"
LIBSTR_RAD_SEQ = "RAD-Seq"
LIBSTR_CLONE = "CLONE"
LIBSTR_POOLCLONE = "POOLCLONE"
LIBSTR_AMPLICON = "AMPLICON"
LIBSTR_CLONEEND = "CLONEEND"
LIBSTR_FINISHING = "FINISHING"
LIBSTR_CHIP_SEQ = "ChIP-Seq"
LIBSTR_MNASE_EQ = "MNase-Seq"
LIBSTR_DNASE_HYPERSENSITIVITY = "DNase-Hypersensitivity"
LIBSTR_BISULFITE_SEQ = "Bisulfite-Seq"
LIBSTR_CTS = "CTS"
LIBSTR_MRE_SEQ = "MRE-Seq"
LIBSTR_MEDIP_SEQ = "MeDIP-Seq"
LIBSTR_MBD_SEQ = "MBD-Seq"
LIBSTR_TN_SEQ = "Tn-Seq"
LIBSTR_VALIDATION = "VALIDATION"
LIBSTR_FAIRE_SEQ = "FAIRE-seq"
LIBSTR_SELEX = "SELEX"
LIBSTR_RIP_SEQ = "RIP-Seq"
LIBSTR_CHIA_PET = "ChIA-PET"
LIBSTR_SYNTHEETIC_LONG_READ = "Synthetic-Long-Read"
LIBSTR_TARGETED_CAPTURE = "Targeted-Capture"
LIBSTR_TETHERED_CHROMATIN_CONFORMATION_CAPTURE = "Tethered Chromatin Conformation Capture"
LIBSTR_OTHER = "OTHER"

#: Allowed library strategies for ENA records.
LIBRARY_STRATEGIES = (
    LIBSTR_WGS,
    LIBSTR_WGA,
    LIBSTR_WXS,
    LIBSTR_RNA_SEQ,
    LIBSTR_SSRNA_SEQ,
    LIBSTR_MIRNA_SEQ,
    LIBSTR_NCRNA_SEQ,
    LIBSTR_FL_CDNA,
    LIBSTR_EST,
    LIBSTR_HIC,
    LIBSTR_ATAC_SEQ,
    LIBSTR_WCS,
    LIBSTR_RAD_SEQ,
    LIBSTR_CLONE,
    LIBSTR_POOLCLONE,
    LIBSTR_AMPLICON,
    LIBSTR_CLONEEND,
    LIBSTR_FINISHING,
    LIBSTR_CHIP_SEQ,
    LIBSTR_MNASE_EQ,
    LIBSTR_DNASE_HYPERSENSITIVITY,
    LIBSTR_BISULFITE_SEQ,
    LIBSTR_CTS,
    LIBSTR_MRE_SEQ,
    LIBSTR_MEDIP_SEQ,
    LIBSTR_MBD_SEQ,
    LIBSTR_TN_SEQ,
    LIBSTR_VALIDATION,
    LIBSTR_FAIRE_SEQ,
    LIBSTR_SELEX,
    LIBSTR_RIP_SEQ,
    LIBSTR_CHIA_PET,
    LIBSTR_SYNTHEETIC_LONG_READ,
    LIBSTR_TARGETED_CAPTURE,
    LIBSTR_TETHERED_CHROMATIN_CONFORMATION_CAPTURE,
    LIBSTR_OTHER,
)

#!/usr/bin/env bash
# Downloads the external tools into ./tools (or the folder given as $1) and prints the tools: block for
# config.yaml. Safe to rerun: tools already in place are skipped. Linux or macOS, needs curl, git, make and g++.
# Tools that need a manual download are listed at the end with instructions.
set -euo pipefail

TOOLS=${1:-tools}
mkdir -p "$TOOLS"
TOOLS=$(cd "$TOOLS" && pwd)
manual=()

# the SimPhy 1.0.2 binaries are x86-64 only
case "$(uname -s)-$(uname -m)" in
    Linux-x86_64) simphy_bin=simphy_lnx64; iqtree_pkg=iqtree-3.0.1-Linux.tar.gz ;;
    Darwin-x86_64) simphy_bin=simphy_mac64; iqtree_pkg=iqtree-3.0.1-macOS.zip ;;
    *) echo "unsupported platform $(uname -s)-$(uname -m), the SimPhy 1.0.2 binaries are x86-64 only" >&2; exit 1 ;;
esac

download() {  # url dest
    if [[ ! -s $2 ]]; then
        echo "downloading $1"
        curl -fL --retry 3 -o "$2.part" "$1"
        mv "$2.part" "$2"
    fi
}

clone() {  # url dest ref
    if [[ ! -d $2/.git ]]; then
        rm -rf "$2.part"
        git clone "$1" "$2.part"
        mv "$2.part" "$2"
    fi
    git -C "$2" checkout --quiet "$3"
}

unpack() {  # archive dest: unpacked into dest.part and renamed, so an interrupted run starts over
    if [[ ! -d $2 ]]; then
        rm -rf "$2.part"
        mkdir -p "$2.part"
        case "$1" in
            *.zip) unzip -q "$1" -d "$2.part" ;;
            *) tar -xzf "$1" -C "$2.part" ;;
        esac
        mv "$2.part" "$2"
    fi
}

need() {  # path that must exist after an install step
    [[ -e $1 ]] || { echo "expected $1 after installing, not found" >&2; exit 1; }
}

# SimPhy 1.0.2
download https://github.com/adamallo/SimPhy/releases/download/v1.0.2/SimPhy_1.0.2.tar.gz "$TOOLS/SimPhy_1.0.2.tar.gz"
unpack "$TOOLS/SimPhy_1.0.2.tar.gz" "$TOOLS/simphy"
simphy=$TOOLS/simphy/SimPhy_1.0.2/bin/$simphy_bin
need "$simphy"
chmod +x "$simphy"

# ASTRAL IV from ASTER v1.23 (source release, built with make)
clone https://github.com/chaoszhang/ASTER.git "$TOOLS/ASTER" v1.23
if [[ ! -x $TOOLS/ASTER/bin/astral4 ]]; then
    make -C "$TOOLS/ASTER" astral
fi
need "$TOOLS/ASTER/bin/astral4"

# IQ-TREE 3.0.1 (AliSim)
download "https://github.com/iqtree/iqtree3/releases/download/v3.0.1/$iqtree_pkg" "$TOOLS/$iqtree_pkg"
unpack "$TOOLS/$iqtree_pkg" "$TOOLS/iqtree-3.0.1"
iqtree=$(find "$TOOLS/iqtree-3.0.1" -type f -name iqtree3 | head -n 1)
[[ -n $iqtree ]] || { echo "no iqtree3 binary under $TOOLS/iqtree-3.0.1" >&2; exit 1; }

# GRANDMA (GRAMPA and GRAMPAIter) at the commit used in the paper, runs on the polybench env python
clone https://github.com/Roshex/Grandma.git "$TOOLS/Grandma" 86abd14
need "$TOOLS/Grandma/gran.py"

# standalone GRAMPA v1.4.1, only for the grampa method's mode: original
clone https://github.com/gwct/grampa.git "$TOOLS/grampa" v1.4.1
need "$TOOLS/grampa/grampa.py"
chmod +x "$TOOLS/grampa/grampa.py"

# Polyphest, upstream main; the commit is printed below so it can be recorded
clone https://github.com/NakhlehLab/Polyphest.git "$TOOLS/Polyphest" main
need "$TOOLS/Polyphest/Polyphest.py"
polyphest_commit=$(git -C "$TOOLS/Polyphest" rev-parse --short HEAD)

# BEAST and TreeAnnotator 1.8.2 (AlloppNET)
download https://github.com/beast-dev/beast-mcmc/releases/download/1.8.2/BEASTv1.8.2.tgz "$TOOLS/BEASTv1.8.2.tgz"
unpack "$TOOLS/BEASTv1.8.2.tgz" "$TOOLS/beast"
beast_dir=$TOOLS/beast/BEASTv1.8.2
need "$beast_dir/bin/beast"
need "$beast_dir/bin/treeannotator"

# MP-SUGAR: PhyNetPy 0.0.8 is installed by environment/mpsugar.yml
if command -v conda > /dev/null && conda env list | grep -q "^mpsugar "; then
    conda run -n mpsugar python -c "import PhyNetPy.MPSugar" \
        || { echo "PhyNetPy.MPSugar does not import in the mpsugar env" >&2; exit 1; }
else
    manual+=("MP-SUGAR: create the env with: conda env create -f environment/mpsugar.yml")
fi

[[ -f $TOOLS/padre-cli.jar ]] || manual+=("PADRE: download the command-line version 2.0 jar from the PADRE site (Lott et al. 2009) and save it as $TOOLS/padre-cli.jar")
[[ -f $TOOLS/PhyloNet_3.8.2.jar ]] || manual+=("PhyloNet (optional, mpallopp_phylonet only): download PhyloNet_3.8.2.jar from the PhyloNet site of the Nakhleh lab and save it as $TOOLS/PhyloNet_3.8.2.jar")
[[ -f $TOOLS/AlloppDT/AlloppDT_5beastxml_toplevel.r && -f $TOOLS/AlloppDT/AlloppDT_6beastxml_bits.r ]] \
    || manual+=("AlloppNET: download the AlloppDT R scripts from Graham Jones's AlloppNET page and put AlloppDT_5beastxml_toplevel.r and AlloppDT_6beastxml_bits.r in $TOOLS/AlloppDT/")

echo
echo "Polyphest commit: $polyphest_commit"
if (( ${#manual[@]} )); then
    echo
    echo "Manual steps:"
    printf '  - %s\n' "${manual[@]}"
fi
echo
echo "Add to config.yaml (tools still to be downloaded are listed too, check their paths):"
cat <<EOF
tools:
  simphy: $simphy
  astral: $TOOLS/ASTER/bin/astral4
  iqtree: $iqtree
  padre_jar: $TOOLS/padre-cli.jar
  grandma: $TOOLS/Grandma/gran.py
  grampa: $TOOLS/grampa/grampa.py
  polyphest: $TOOLS/Polyphest/Polyphest.py
  phylonet_jar: $TOOLS/PhyloNet_3.8.2.jar
  beast: $beast_dir/bin/beast
  treeannotator: $beast_dir/bin/treeannotator
  alloppdt: $TOOLS/AlloppDT
EOF

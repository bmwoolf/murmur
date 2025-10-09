#!/bin/bash

# detect system information and generates a personalized config.yml

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
TEMPLATE_FILE="$SCRIPT_DIR/config.yml.example"
OUTPUT_FILE="$PROJECT_ROOT/config.yml"

echo "~~~detecting system information~~~"

# OS information
OS_NAME=$(lsb_release -si 2>/dev/null || echo "Linux")
OS_VERSION=$(lsb_release -sr 2>/dev/null || echo "Unknown")
KERNEL_VERSION=$(uname -r)
ARCHITECTURE=$(uname -m)

# CPU information
CPU_MODEL=$(lscpu | grep "Model name" | cut -d: -f2 | xargs)
PHYSICAL_CORES=$(lscpu | grep "^CPU(s):" | cut -d: -f2 | xargs)
LOGICAL_CORES=$(lscpu | grep "^CPU(s):" | cut -d: -f2 | xargs)
CPU_FREQUENCY=$(lscpu | grep "CPU MHz" | cut -d: -f2 | xargs | cut -d. -f1)

# memory information (in GB)
TOTAL_RAM_GB=$(free -g | grep "Mem:" | awk '{print $2}')
AVAILABLE_RAM_GB=$(free -g | grep "Mem:" | awk '{print $7}')

# storage information (in GB)
TOTAL_STORAGE_GB=$(df -BG / | tail -1 | awk '{print $2}' | sed 's/G//')
FREE_STORAGE_GB=$(df -BG / | tail -1 | awk '{print $4}' | sed 's/G//')

# GPU information
if command -v nvidia-smi &> /dev/null; then
    GPU_MODEL=$(nvidia-smi --query-gpu=name --format=csv,noheader,nounits | head -1)
    GPU_VRAM_GB=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1 | awk '{print int($1/1024)}')
    GPU_COMPUTE_CAPABILITY=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader,nounits | head -1)
    GPU_MULTIPROCESSORS=$(nvidia-smi --query-gpu=multiprocessor_count --format=csv,noheader,nounits | head -1)
    CUDA_VERSION=$(nvidia-smi --query-gpu=driver_version --format=csv,noheader,nounits | head -1)
else
    GPU_MODEL="None detected"
    GPU_VRAM_GB="0"
    GPU_COMPUTE_CAPABILITY="N/A"
    GPU_MULTIPROCESSORS="0"
    CUDA_VERSION="N/A"
fi

# Python/PyTorch information
if command -v python3 &> /dev/null; then
    PYTHON_VERSION=$(python3 --version | cut -d' ' -f2)
else
    PYTHON_VERSION="Not installed"
fi

PYTORCH_VERSION="Not installed"
if command -v python3 &> /dev/null; then
    PYTORCH_VERSION=$(python3 -c "import torch; print(torch.__version__)" 2>/dev/null || echo "Not installed")
fi

# pipeline configuration (calculated based on system specs)
ANNOTATION_THREADS=$((LOGICAL_CORES > 4 ? LOGICAL_CORES - 4 : LOGICAL_CORES))
REFERENCE_GENOME="GRCh38"
VEP_CACHE_DIR="~/.vep"
SNPEFF_DATA_DIR="~/.snpeff/data"
VEP_DOCKER_IMAGE="ensemblorg/ensembl-vep:release_109.3"
SNPEFF_DOCKER_IMAGE="staphb/snpeff:latest"

# X-Atlas streaming configuration
XATLAS_H5AD="null"  # null = streaming mode (recommended)
XATLAS_CELL_LINE="HCT116"  # HCT116 or HEK293T

# CPA model settings based on available resources
if [ "$GPU_VRAM_GB" -gt 8 ]; then
    CPA_BATCH_SIZE="512"
elif [ "$GPU_VRAM_GB" -gt 4 ]; then
    CPA_BATCH_SIZE="256"
else
    CPA_BATCH_SIZE="128"
fi

CPA_LEARNING_RATE="0.001"
CPA_MAX_EPOCHS="100"

# Memory Management (leave headroom)
MAX_RAM_USAGE_GB=$((AVAILABLE_RAM_GB > 4 ? AVAILABLE_RAM_GB - 4 : AVAILABLE_RAM_GB))
GPU_MEMORY_FRACTION="0.9"

# Output Settings
KEEP_INTERMEDIATE_FILES="true"
KEEP_TEMP_FILES="false"
COMPRESSION_TYPE="gzip"

# Check if template exists
if [ ! -f "$TEMPLATE_FILE" ]; then
    echo "Error: Template file $TEMPLATE_FILE not found!"
    exit 1
fi

# Create backup of existing config if it exists
if [ -f "$OUTPUT_FILE" ]; then
    echo "~~~creating backup of existing config.yml~~~"
    cp "$OUTPUT_FILE" "$OUTPUT_FILE.backup.$(date +%Y%m%d_%H%M%S)"
fi

# Generate the config file by replacing placeholders
sed -e "s/{{OS_NAME}}/$OS_NAME/g" \
    -e "s/{{OS_VERSION}}/$OS_VERSION/g" \
    -e "s/{{KERNEL_VERSION}}/$KERNEL_VERSION/g" \
    -e "s/{{ARCHITECTURE}}/$ARCHITECTURE/g" \
    -e "s/{{CPU_MODEL}}/$CPU_MODEL/g" \
    -e "s/{{PHYSICAL_CORES}}/$PHYSICAL_CORES/g" \
    -e "s/{{LOGICAL_CORES}}/$LOGICAL_CORES/g" \
    -e "s/{{CPU_FREQUENCY}}/$CPU_FREQUENCY/g" \
    -e "s/{{TOTAL_RAM_GB}}/$TOTAL_RAM_GB/g" \
    -e "s/{{AVAILABLE_RAM_GB}}/$AVAILABLE_RAM_GB/g" \
    -e "s/{{TOTAL_STORAGE_GB}}/$TOTAL_STORAGE_GB/g" \
    -e "s/{{FREE_STORAGE_GB}}/$FREE_STORAGE_GB/g" \
    -e "s/{{GPU_MODEL}}/$GPU_MODEL/g" \
    -e "s/{{GPU_VRAM_GB}}/$GPU_VRAM_GB/g" \
    -e "s/{{GPU_COMPUTE_CAPABILITY}}/$GPU_COMPUTE_CAPABILITY/g" \
    -e "s/{{GPU_MULTIPROCESSORS}}/$GPU_MULTIPROCESSORS/g" \
    -e "s/{{CUDA_VERSION}}/$CUDA_VERSION/g" \
    -e "s/{{PYTHON_VERSION}}/$PYTHON_VERSION/g" \
    -e "s/{{PYTORCH_VERSION}}/$PYTORCH_VERSION/g" \
    -e "s/{{ANNOTATION_THREADS}}/$ANNOTATION_THREADS/g" \
    -e "s/{{REFERENCE_GENOME}}/$REFERENCE_GENOME/g" \
    -e "s|{{VEP_CACHE_DIR}}|$VEP_CACHE_DIR|g" \
    -e "s|{{SNPEFF_DATA_DIR}}|$SNPEFF_DATA_DIR|g" \
    -e "s|{{VEP_DOCKER_IMAGE}}|$VEP_DOCKER_IMAGE|g" \
    -e "s|{{SNPEFF_DOCKER_IMAGE}}|$SNPEFF_DOCKER_IMAGE|g" \
    -e "s/{{CPA_BATCH_SIZE}}/$CPA_BATCH_SIZE/g" \
    -e "s/{{CPA_LEARNING_RATE}}/$CPA_LEARNING_RATE/g" \
    -e "s/{{CPA_MAX_EPOCHS}}/$CPA_MAX_EPOCHS/g" \
    -e "s/{{MAX_RAM_USAGE_GB}}/$MAX_RAM_USAGE_GB/g" \
    -e "s/{{GPU_MEMORY_FRACTION}}/$GPU_MEMORY_FRACTION/g" \
    -e "s/{{KEEP_INTERMEDIATE_FILES}}/$KEEP_INTERMEDIATE_FILES/g" \
    -e "s/{{KEEP_TEMP_FILES}}/$KEEP_TEMP_FILES/g" \
    -e "s/{{COMPRESSION_TYPE}}/$COMPRESSION_TYPE/g" \
    "$TEMPLATE_FILE" > "$OUTPUT_FILE"

echo "Successfully generated config.yml at $OUTPUT_FILE"
echo ""
echo "System Summary:"
echo "   OS: $OS_NAME $OS_VERSION"
echo "   CPU: $CPU_MODEL ($LOGICAL_CORES cores)"
echo "   RAM: ${TOTAL_RAM_GB}GB total, ${AVAILABLE_RAM_GB}GB available"
echo "   GPU: $GPU_MODEL (${GPU_VRAM_GB}GB VRAM)"
echo "   Python: $PYTHON_VERSION"
echo "   PyTorch: $PYTORCH_VERSION"
echo ""
echo "Pipeline Settings:"
echo "   Annotation threads: $ANNOTATION_THREADS"
echo "   CPA batch size: $CPA_BATCH_SIZE"
echo "   Max RAM usage: ${MAX_RAM_USAGE_GB}GB"
echo ""
echo "You can now edit $OUTPUT_FILE to customize settings as needed."

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdint>
#include <limits>
#include <algorithm>

const double PI = 3.14159265358979323846;
const int MAX_OFFSET = 5000; // Limit the offset search range
const double COMPRESS_THRESHOLD = 30000.0; // Dynamic range compression threshold
const double COMPRESS_RATIO = 2.0; // Compression ratio

enum BlendMode {
    AVERAGE,
    CROSSFADE,
    ADDITIVE,
    SUBTRACTIVE
};

// WAV file header
struct WAVHeader {
    char riff[4];
    uint32_t fileSize;
    char wave[4];
    char fmt[4];
    uint32_t fmtSize;
    uint16_t audioFormat;
    uint16_t numChannels;
    uint32_t sampleRate;
    uint32_t byteRate;
    uint16_t blockAlign;
    uint16_t bitsPerSample;
    char data[4];
    uint32_t dataSize;
};

// Read WAV file (supporting stereo)
bool readWAV(const std::string& filename, WAVHeader& header, std::vector<int16_t>& samples) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) return false;

    file.read(reinterpret_cast<char*>(&header), sizeof(WAVHeader));

    // Calculate number of samples for stereo (2 channels per sample)
    size_t numSamples = header.dataSize / (header.bitsPerSample / 8) / header.numChannels;
    samples.resize(numSamples * header.numChannels);  // Stereo: left and right channel per sample

    file.read(reinterpret_cast<char*>(samples.data()), header.dataSize);
    file.close();
    return true;
}

// Write WAV file (supporting stereo)
bool writeWAV(const std::string& filename, WAVHeader header, const std::vector<int16_t>& samples) {
    header.dataSize = samples.size() * sizeof(int16_t);
    header.fileSize = header.dataSize + sizeof(WAVHeader) - 8;

    std::ofstream file(filename, std::ios::binary);
    if (!file) return false;

    file.write(reinterpret_cast<const char*>(&header), sizeof(WAVHeader));
    file.write(reinterpret_cast<const char*>(samples.data()), header.dataSize);
    file.close();
    return true;
}

// Align and blend two audio signals (Stereo)
std::vector<int16_t> alignAudio(const std::vector<int16_t>& samples1, const std::vector<int16_t>& samples2, BlendMode mode, uint16_t numChannels, double crossfadePercentage) {
    size_t numSamples = std::min(samples1.size(), samples2.size()) / numChannels;
    std::vector<int16_t> aligned(numSamples * numChannels);

    for (size_t i = 0; i < numSamples; ++i) {
        for (uint16_t ch = 0; ch < numChannels; ++ch) {
            size_t index = i * numChannels + ch;
            switch (mode) {
            case AVERAGE:
                aligned[index] = (samples1[index] + samples2[index]) / 2;
                break;
            case CROSSFADE:
            {
                // Adjust the blend factor based on crossfadePercentage
                double blendFactor = static_cast<double>(i) / numSamples;
                // Apply crossfade adjustment (user-defined percentage of fade)
                blendFactor *= crossfadePercentage;
                aligned[index] = static_cast<int16_t>((samples1[index] * (1.0 - blendFactor)) + (samples2[index] * blendFactor));
            }
            break;
            case ADDITIVE:
                aligned[index] = std::max(-32768, std::min(32767, samples1[index] + samples2[index]));
                break;
            case SUBTRACTIVE:
                aligned[index] = std::max(-32768, std::min(32767, samples1[index] - samples2[index]));
                break;
            }
        }
    }
    return aligned;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input1.wav> <input2.wav> <output.wav>" << std::endl;
        return 1;
    }
    std::string file1 = argv[1];
    std::string file2 = argv[2];
    std::string outputFile = argv[3];

    WAVHeader header1, header2;
    std::vector<int16_t> samples1, samples2;

    if (!readWAV(file1, header1, samples1) || !readWAV(file2, header2, samples2)) {
        std::cerr << "Error: Could not read input files." << std::endl;
        return 1;
    }

    // Check if both files have the same format (sample rate, bit depth, and number of channels)
    if (header1.sampleRate != header2.sampleRate || header1.bitsPerSample != header2.bitsPerSample ||
        header1.numChannels != header2.numChannels) {
        std::cerr << "Error: WAV files must have the same format." << std::endl;
        return 1;
    }

    int modeSelection;
    std::cout << "Select blending mode:\n1. Average\n2. Crossfade\n3. Additive\n4. Subtractive\nChoice: ";
    std::cin >> modeSelection;
    BlendMode mode = static_cast<BlendMode>(modeSelection - 1);

    double crossfadePercentage = 1.0;  // Default 100% crossfade

    // If Crossfade is selected, allow the user to adjust the blend percentage
    if (mode == CROSSFADE) {
        std::cout << "Enter crossfade percentage (0.0 to 1.0, default 1.0): ";
        std::cin >> crossfadePercentage;
        crossfadePercentage = std::max(0.0, std::min(1.0, crossfadePercentage));  // Ensure the percentage is between 0 and 1
    }

    std::vector<int16_t> alignedSamples = alignAudio(samples1, samples2, mode, header1.numChannels, crossfadePercentage);

    if (!writeWAV(outputFile, header1, alignedSamples)) {
        std::cerr << "Error: Could not write output file." << std::endl;
        return 1;
    }

    std::cout << "Aligned WAV saved as: " << outputFile << std::endl;
    return 0;
}

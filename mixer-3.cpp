#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <complex>

const double PI = 3.14159265358979323846;
using Complex = std::complex<double>;

// WAV Header structure (assuming 16-bit PCM)
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

// Function to read a WAV file and extract left and right channels
bool readWAV(const std::string& filename, std::vector<double>& left, std::vector<double>& right, WAVHeader& header) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return false;
    }

    file.read(reinterpret_cast<char*>(&header), sizeof(WAVHeader));
    if (header.audioFormat != 1 || header.bitsPerSample != 16) {
        std::cerr << "Unsupported WAV format! Only 16-bit PCM is supported." << std::endl;
        return false;
    }

    int numSamples = header.dataSize / (header.bitsPerSample / 8);
    left.resize(numSamples / header.numChannels);
    right.resize(numSamples / header.numChannels);

    for (int i = 0; i < numSamples / header.numChannels; ++i) {
        int16_t sampleL, sampleR;
        file.read(reinterpret_cast<char*>(&sampleL), sizeof(sampleL));
        left[i] = sampleL / 32768.0; // Normalize to -1.0 to 1.0

        if (header.numChannels == 2) {
            file.read(reinterpret_cast<char*>(&sampleR), sizeof(sampleR));
            right[i] = sampleR / 32768.0;
        }
    }

    file.close();
    return true;
}

// Function to write a stereo WAV file
bool writeWAV(const std::string& filename, const std::vector<double>& left, const std::vector<double>& right, const WAVHeader& header) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Could not open " << filename << " for writing." << std::endl;
        return false;
    }

    WAVHeader outHeader = header;
    outHeader.dataSize = left.size() * 2 * header.numChannels;
    outHeader.fileSize = outHeader.dataSize + sizeof(WAVHeader) - 8;

    file.write(reinterpret_cast<const char*>(&outHeader), sizeof(WAVHeader));

    for (size_t i = 0; i < left.size(); ++i) {
        int16_t sampleL = static_cast<int16_t>(std::max(-1.0, std::min(1.0, left[i])) * 32767);
        file.write(reinterpret_cast<const char*>(&sampleL), sizeof(sampleL));

        if (header.numChannels == 2) {
            int16_t sampleR = static_cast<int16_t>(std::max(-1.0, std::min(1.0, right[i])) * 32767);
            file.write(reinterpret_cast<const char*>(&sampleR), sizeof(sampleR));
        }
    }

    file.close();
    return true;
}

// FFT functions
void fft(std::vector<Complex>& x) {
    int N = x.size();
    if (N <= 1) return;

    std::vector<Complex> even(N / 2), odd(N / 2);
    for (int i = 0; i < N / 2; ++i) {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }

    fft(even);
    fft(odd);

    for (int k = 0; k < N / 2; ++k) {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k + N / 2] = even[k] - t;
    }
}

void ifft(std::vector<Complex>& X) {
    int N = X.size();
    for (auto& x : X) x = std::conj(x);
    fft(X);
    for (auto& x : X) x = std::conj(x) / (double)N;
}

// Hanning Window
void applyHanningWindow(std::vector<double>& frame) {
    int N = frame.size();
    for (int i = 0; i < N; ++i) {
        frame[i] *= 0.5 * (1 - std::cos(2 * PI * i / (N - 1)));
    }
}

// Spectral mixing with mixing ratio
void spectralMix(const std::vector<Complex>& X1, const std::vector<Complex>& X2, std::vector<Complex>& Y, double mixRatio) {
    for (size_t i = 0; i < X1.size(); ++i) {
        double mag1 = std::abs(X1[i]);
        double mag2 = std::abs(X2[i]);
        double phase = std::arg(X1[i]); // Preserve phase of first signal

        // Apply the mixing ratio to the magnitudes
        double newMag = (mixRatio * mag1) + ((1.0 - mixRatio) * mag2);
        Y[i] = std::polar(newMag, phase);
    }
}

// Process spectral mixing for stereo, handle both left and right channels separately
std::vector<double> processSpectralMixing(const std::vector<double>& left1, const std::vector<double>& right1,
    const std::vector<double>& left2, const std::vector<double>& right2,
    int frameSize, int hopSize, double mixRatio) {
    size_t numFrames = (left1.size() - frameSize) / hopSize;
    std::vector<double> mixedLeft(left1.size(), 0.0);
    std::vector<double> mixedRight(right1.size(), 0.0);

    for (size_t i = 0; i < numFrames; ++i) {
        std::vector<double> frameLeft1(left1.begin() + i * hopSize, left1.begin() + i * hopSize + frameSize);
        std::vector<double> frameRight1(right1.begin() + i * hopSize, right1.begin() + i * hopSize + frameSize);
        std::vector<double> frameLeft2(left2.begin() + i * hopSize, left2.begin() + i * hopSize + frameSize);
        std::vector<double> frameRight2(right2.begin() + i * hopSize, right2.begin() + i * hopSize + frameSize);

        // Apply Hanning window to both channels
        applyHanningWindow(frameLeft1);
        applyHanningWindow(frameRight1);
        applyHanningWindow(frameLeft2);
        applyHanningWindow(frameRight2);

        // Perform FFT on both frames for both channels
        std::vector<Complex> X1Left(frameLeft1.begin(), frameLeft1.end());
        std::vector<Complex> X1Right(frameRight1.begin(), frameRight1.end());
        std::vector<Complex> X2Left(frameLeft2.begin(), frameLeft2.end());
        std::vector<Complex> X2Right(frameRight2.begin(), frameRight2.end());

        fft(X1Left);
        fft(X1Right);
        fft(X2Left);
        fft(X2Right);

        // Mix the spectra for both channels with the mixing ratio
        std::vector<Complex> YLeft(frameSize), YRight(frameSize);
        spectralMix(X1Left, X2Left, YLeft, mixRatio);
        spectralMix(X1Right, X2Right, YRight, mixRatio);

        // Perform IFFT to get time-domain signals back
        ifft(YLeft);
        ifft(YRight);

        // Overlap-add into the output signals
        for (size_t j = 0; j < frameSize; ++j) {
            mixedLeft[i * hopSize + j] += YLeft[j].real();
            mixedRight[i * hopSize + j] += YRight[j].real();
        }
    }

    return mixedLeft; // You would similarly return mixedRight in your main function
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " input1.wav input2.wav output.wav" << std::endl;
        return 1;
    }

    std::string inputFile1 = argv[1];
    std::string inputFile2 = argv[2];
    std::string outputFile = argv[3];

    WAVHeader header1, header2;
    std::vector<double> left1, right1, left2, right2;

    if (!readWAV(inputFile1, left1, right1, header1) || !readWAV(inputFile2, left2, right2, header2)) {
        return 1;
    }

    std::cout << "Performing Spectral Mixing (Slow)." << std::endl << std::endl;

    // Prompt for mixing ratio
    double mixRatio;
    std::cout << "Enter mixing ratio (0.0 to 1.0): ";
    std::cin >> mixRatio;

    std::vector<double> mixedLeft = processSpectralMixing(left1, right1, left2, right2, 1024, 512, mixRatio);
    std::vector<double> mixedRight = processSpectralMixing(left1, right1, left2, right2, 1024, 512, mixRatio);

    if (!writeWAV(outputFile, mixedLeft, mixedRight, header1)) {
        return 1;
    }

    std::cout << "Stereo spectral mixing complete! Output saved to " << outputFile << std::endl;
    return 0;
}

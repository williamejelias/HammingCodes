# Hamming Encoder-Decoder Simulator

This code simulates the encoding, transmission over a noisy channel, and decoding of a bit string.

## Usage

Inputs:
* message: Array of binary values e.g. [0,1,1,0]
* bsc: Binary Symmetric Channel (BSC) Probability e.g. 0.5

The BSC value denotes the probability that each bit is flipped in transmission in order to simulate a noisy channel.

```bash
python3 hammingSimulator.py message bsc
```

The file outputs the values decoded through three decoding techniques:
* Brute Force
    * Search all valid codewords and their hamming distance from the received vector - decode the valid codeword with the smallest distance.
* Local Search
    * Search for the valid codeword in the Hamming neighbourhood by flipping a bit at an index and attempting to decode.
* Syndrome
    * Syndrome decoding is minimum distance decoding using a reduced lookup table. This is allowed by the linearity of the code.

The output shows whether decoding was successful or not, and whether or not the three decoding algorithms agreed with each others decoding.

Refer to: 
https://en.wikipedia.org/wiki/Hamming_code
https://en.wikipedia.org/wiki/Decoding_methods
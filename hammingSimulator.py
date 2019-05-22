import numpy as np
import random
import itertools
import sys


# function HammingG
# input: a number r
# output: G, the generator matrix of the (2^r-1,2^r-r-1) Hamming code
def hamming_generator_matrix(r):
    n = 2 ** r - 1

    # construct permutation pi
    pi = []
    for i in range(r):
        pi.append(2 ** (r - i - 1))
    for j in range(1, r):
        for k in range(2 ** j + 1, 2 ** (j + 1)):
            pi.append(k)

    # construct rho = pi^(-1)
    rho = []
    for i in range(n):
        rho.append(pi.index(i + 1))

    # construct H'
    H = []
    for i in range(r, n):
        H.append(decimal_to_vector(pi[i], r))

    # construct G'
    GG = [list(i) for i in zip(*H)]
    for i in range(n - r):
        GG.append(decimal_to_vector(2 ** (n - r - i - 1), n - r))

    # apply rho to get Gtranpose
    G = []
    for i in range(n):
        G.append(GG[rho[i]])

    # transpose
    G = [list(i) for i in zip(*G)]

    return G


# function decimalToVector
# input: numbers n and r (0 <= n<2**r)
# output: a string v of r bits representing n
def decimal_to_vector(n, r):
    v = []
    for s in range(r):
        v.insert(0, n % 2)
        n //= 2
    return v


# function message
# input: vector a
# output: the encoded mesage from the vector
def message(a):
    v = []
    r = 2
    while (2 ** r) - 2 * r - 1 < len(a):
        r += 1
    k = ((2 ** r) - r - 1)
    for i in decimal_to_vector(len(a), r):
        v.append(i)
    for j in a:
        v.append(j)
    while len(v) < k:
        v.append(0)
    return v


# function hammingEncoder
# input: message m
# output: hamming code of message (error if m is not of length (2**r - r - 1) for some r >= 2)
def hamming_encoder(m):
    h = []
    # error check for length
    r = 2
    valid_length = (2 ** r) - r - 1
    while valid_length < len(m):
        r += 1
        valid_length = ((2 ** r) - r - 1)
    if valid_length > len(m):
        print("r: ", r)
        string = "Error: incorrect length"
        return string
    elif valid_length == len(m):
        # rest of function
        a = np.array(m)
        b = np.array(hamming_generator_matrix(r))
        c = np.dot(a, b)
        for i in c:
            if i > 1:
                if i % 2 == 0:
                    h.append(0)
                else:
                    h.append(1)
            else:
                h.append(i)
        return h


# function BSC
# input: vector of any length c, probability p
# output: code with probability of bit flipping p
def BSC(c, p):
    v = []
    for i in c:
        x = random.randint(1, 100) / 100
        if x < p:
            if i == 1:
                v.append(0)
            else:
                v.append(1)
        else:
            v.append(i)
    return v


# function hammingBruteForce
# input: vector v where len(v) is (2**r - 1) for some r >= 2
# output: brute force hamming code
def hamming_brute_force(v):
    actual_codeword = ''
    r = 2
    valid_length = 2 ** r - 1
    while 2 ** r - 1 < len(v):
        r += 1
        valid_length = 2 ** r - 1
    if valid_length > len(v):
        print("r: ", r)
        string = "Error: incorrect length"
        return string
    elif valid_length == len(v):
        # rest of function
        # generate all bit_strings of length 2**r - r - 1
        bit_strings = [''.join(i) for i in itertools.product('01', repeat=(2 ** r - r - 1))]
        list_of_bitstrings = []
        for i in bit_strings:
            bitstring = []
            for char in i:
                bitstring.append(int(char))
            list_of_bitstrings.append(bitstring)
        # generate the codeword associated with each bitstring and store in a list
        list_of_generated_codewords = []
        b = np.array(hamming_generator_matrix(r))
        for i in list_of_bitstrings:
            a = np.array(i)
            c = np.dot(a, b)
            cw = []
            for i in c:
                if i > 1:
                    if i % 2 == 0:
                        cw.append(0)
                    else:
                        cw.append(1)
                else:
                    cw.append(i)
            list_of_generated_codewords.append(cw)
        # find the codeword in the list which has the smallest distance from the input vector
        minimum_difference = len(v)
        print("Decoding by brute force")
        for i in range(len(list_of_generated_codewords)):
            currentDifference = 0
            for bitIndex in range(len(v)):
                if v[bitIndex] != list_of_generated_codewords[i][bitIndex]:
                    currentDifference += 1
            if currentDifference < minimum_difference:
                minimum_difference = currentDifference
                actual_codeword = list_of_generated_codewords[i]
            # print("m", i, " G = ", list_of_generated_codewords[i])
            # print("d(v, m", i, " G) = ", currentDifference)
            if minimum_difference == 1:
                break
            # else:
            #     print()
        print("hatc: ", actual_codeword, "  -  Minimum Difference: ", minimum_difference)
    return actual_codeword


# function hammingLocalSearch
# input: vector v where len(v) is (2**r - 1) for some r >= 2
# output: local search hamming code
def hamming_local_search(v):
    r = 2
    valid_length = 2 ** r - 1
    while 2 ** r - 1 < len(v):
        r += 1
        valid_length = 2 ** r - 1
    if valid_length > len(v):
        print("r: ", r)
        string = "Error: incorrect length"
        return string
    elif valid_length == len(v):
        print("Decoding by local search")
        # rest of function
        # build parity check matrix, which has dimensions r rows * (2**r - 1) columns
        # e.g [[0, 1, 1], [1, 0, 1]], or [[0, 0, 0, 1, 1, 1, 1], [0, 1, 1, 0, 0, 1, 1], [1, 0, 1, 0, 1, 0, 1]]
        list_of_binary_numbers = []
        for i in range(2 ** r - 1):
            list_of_binary_numbers.append(decimal_to_vector(i + 1, r))
        # row 0 comprises of the element at index 0 of the first 2**r - 1 binary digits
        parity_check_matrix = []
        for j in range(r):
            row = []
            for k in list_of_binary_numbers:
                row.append(k[j])
            parity_check_matrix.append(row)
        # if v is a codeword, decode using parity check matrix
        a = np.array(v)
        b = np.array(parity_check_matrix)
        c = b.T
        d = np.dot(a, c)
        decoded = []
        for l in d:
            if l > 1:
                if l % 2 == 0:
                    decoded.append(0)
                else:
                    decoded.append(1)
            else:
                decoded.append(l)
        all_zeroes = not np.array(decoded).any()
        if all_zeroes:
            print("syndrome = ", decoded)
            print("hatc: ", v)
            return v
        else:
            # else, v with some bit switch will be a codeword (d_min is 1)
            list_of_iterations_of_input = []
            for m in range(0, len(v)):
                vCopy = list(v)
                if vCopy[m] == 0:
                    vCopy[m] = 1
                elif vCopy[m] == 1:
                    vCopy[m] = 0
                list_of_iterations_of_input.append(vCopy)
            list_of_decoded_inputs = []
            for n in range(len(list_of_iterations_of_input)):
                e = np.array(list_of_iterations_of_input[n])
                f = np.dot(e, c)
                decoded_iteration = []
                for j in f:
                    if j > 1:
                        if j % 2 == 0:
                            decoded_iteration.append(0)
                        else:
                            decoded_iteration.append(1)
                    else:
                        decoded_iteration.append(j)
                list_of_decoded_inputs.append(decoded_iteration)
                print("v + e", n + 1, " = ", list_of_iterations_of_input[n])
                print("syndrome = ", decoded_iteration)
                if not np.array(decoded_iteration).any():
                    print("hatc = ", list_of_iterations_of_input[n])
                    return list_of_iterations_of_input[n]
                else:
                    print()


# function hammingSyndrome
# input: vector v where len(v) is (2**r - 1) for some r >= 2
# output: syndrome search hamming code
def hamming_syndrome(v):
    r = 2
    valid_length = 2 ** r - 1
    while 2 ** r - 1 < len(v):
        r += 1
        valid_length = 2 ** r - 1
    if valid_length > len(v):
        print("r: ", r)
        string = "Error: incorrect length"
        return string
    elif valid_length == len(v):
        print("Decoding by syndrome")
        # rest of function
        # build parity check matrix, which has dimensions r rows * (2**r - 1) columns
        # e.g [[0, 1, 1], [1, 0, 1]], or [[0, 0, 0, 1, 1, 1, 1], [0, 1, 1, 0, 0, 1, 1], [1, 0, 1, 0, 1, 0, 1]]
        list_of_binary_numbers = []
        for i in range(2 ** r - 1):
            list_of_binary_numbers.append(decimal_to_vector(i + 1, r))
        # row 0 comprises of the element at index 0 of the first 2**r - 1 binary digits
        parity_check_matrix = []
        for j in range(r):
            row = []
            for k in list_of_binary_numbers:
                row.append(k[j])
            parity_check_matrix.append(row)
        # if v is a codeword, decode using parity check matrix
        a = np.array(v)
        b = np.array(parity_check_matrix)
        c = b.T
        d = np.dot(a, c)
        decoded = []
        for l in d:
            if l > 1:
                if l % 2 == 0:
                    decoded.append(0)
                else:
                    decoded.append(1)
            else:
                decoded.append(l)
        decoded_vector = np.array(decoded)
        all_zeroes = not decoded_vector.any()
        if all_zeroes:
            print("syndrome = ", decoded)
            print("hatc = ", v)
            return v
        else:
            decoded_string = ''
            for i in decoded:
                decoded_string += str(i)
            decoded_reversed = ''.join(reversed(decoded_string))
            value = 0
            for i in range(len(decoded_reversed)):
                value += 2 ** i * int(decoded_reversed[i])
            if v[value - 1] == 1:
                v[value - 1] = 0
            elif v[value - 1] == 0:
                v[value - 1] = 1
            print("syndrome = ", decoded)
            print("i = ", value)
            print("hatc = ", v)
            return v


# function messageFromCodeword
# input: codeword c where len(c) is (2**r - 1) for some r >= 2
# output: message
def message_from_codeword(c):
    r = 2
    valid_length = 2 ** r - 1
    while 2 ** r - 1 < len(c):
        r += 1
        valid_length = 2 ** r - 1
    if valid_length > len(c):
        string = "Error: incorrect length"
        return string
    elif valid_length == len(c):
        # rest of function
        word_to_message = c
        number_of_values_removed = 0
        for i in range(r):
            del word_to_message[2 ** i - 1 - number_of_values_removed]
            number_of_values_removed += 1
        return word_to_message


# function dataFromMessage
# input: binary message m where len(m) is (2**r - r - 1) for some r >= 2
# output: raw data
def data_from_message(m):
    data = []
    r = 2
    valid_length = (2 ** r) - r - 1
    while valid_length < len(m):
        r += 1
        valid_length = ((2 ** r) - r - 1)
    if valid_length > len(m):
        string = "Error: incorrect length"
        return string
    elif valid_length == len(m):
        # rest of function
        l_binary = []
        for i in range(r):
            l_binary.append(m[0])
            del m[0]
        l_binary_string = ''
        for i in l_binary:
            l_binary_string += str(i)
        value = int(l_binary_string, 2)
        for i in range(value):
            if i > len(m) - 1:
                return data
            else:
                data.append(m[i])
        return data


# function simulation
# input: raw data a, probability p, where 0 <= p <= 1
# output: string denoting whether the transmission was successful or not
def simulation(a, p):
    print("Raw data")
    print("a: ", a)
    m = message(a)
    print("\nMessage of raw data")
    print("m: ", m)
    hEncoded = hamming_encoder(m)
    print("\nHamming Encoded Codeword")
    print("c: ", hEncoded)
    noisified = BSC(hEncoded, p)
    print("\nReceived vector after noisy transmission")
    print("v: ", noisified)
    print()
    decoded1 = hamming_brute_force(noisified)
    print()
    decoded2 = hamming_local_search(noisified)
    print()
    decoded3 = hamming_syndrome(noisified)
    if decoded1 == decoded2 == decoded3:
        print()
        print("Deduced Codeword")
        print("hatc: ", decoded1)
        print()
        decodedMessage = message_from_codeword(decoded1)
        print("Converted to Message")
        print("hatm: ", decodedMessage)
        print()
        decodedRawData = data_from_message(decodedMessage)
        print("Received Raw data")
        print("hata: ", decodedRawData)
        if a == decodedRawData:
            print()
            print("Success! - All decoding methods return the same data.")
        else:
            print()
            print("Transmission failure - Received decoded data is different to transmitted data ")
    else:
        print("Decoder error - decoded message is different from different decoding algorithms")


# get command line data
try:
    r_data = sys.argv[1]
    bsc = sys.argv[2]
    # check that command line data is of the right format
    try:
        raw_data = []
        for i in r_data:
            if i == "1":
                raw_data.append(1)
            elif i == "0":
                raw_data.append(0)
        print(r_data, raw_data, type(r_data))
        simulation(raw_data, float(bsc))
    except Exception as msg:
        print(msg)
except Exception:
    print("Not enough arguments")
    exit()



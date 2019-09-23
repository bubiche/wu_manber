# Wu - Manber Implementation (C++)
This is a header only implementation of the Wu-Manber algorithm for multiple pattern matching: http://webglimpse.net/pubs/TR94-17.pdf, I cannot find any C++ implementation on github so I decided to make one.

This algorithm is one of the ones that power agrep, the use of hash table is quite interesting since most popular algorithms use some sort of finite automata structure.

My implementation use some C++11 features, I tried my best to make it not a mess of C mixed with C++
Usage example can be found in main.cpp

```c++
wu_manber::WuManber<char> wuManber;
std::vector<std::string> patternList = { "some", "word", "abracadabra" };
std::string text = "someabracadabra";
wuManber.preProcess(patternList);
wuManber.scan(text,
  [](const std::string &pattern, size_t indexInPatternList, size_t startIndexInText) {
    std::cout << "Matched pattern #" << indexInPatternList <<": '" << pattern << "' at position " << startIndexInText << std::endl;
  }
);
```

Currently only works with patterns of length >= 3. For patterns of length 1 and 2, we could use 2 arrays `lengthOneLookup[ALPHABET_SIZE]` and `lengthTwoLookup[ALPHABET_SIZE][ALPHABET_SIZE]` to scan for them as described in this paper: https://www.researchgate.net/publication/275266493_Improve_Wu-Manber_Algorithm_for_Multiple_Pattern_Matching_Time_Efficiency_Optimization

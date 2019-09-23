#include "wu_manber.hpp"
#include <string>
#include <vector>
#include <iostream>

int main()
{
  // simple driver program for th Wu - Manber algorithm
  wu_manber::WuManber<char> wuManber;
  std::vector<std::string> patternList = { "some", "word", "abracadabra", "x", "ra" };
  std::string text = "xrayx someabracadabra";
  wuManber.preProcess(patternList);
  wuManber.scan(text,
    [](const std::string &pattern, size_t indexInPatternList, size_t startIndexInText) {
      std::cout << "Matched pattern #" << indexInPatternList <<": '" << pattern << "' at position " << startIndexInText << std::endl;
    }
  );
  return 0;
}

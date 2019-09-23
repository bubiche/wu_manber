#ifndef WU_MANBER_HPP
#define WU_MANBER_HPP

#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <functional>

namespace wu_manber {

  namespace { // anonymous namespace, things in here are "private" to wu_manber namespace

    // fast mod (ref: https://www.youtube.com/watch?v=nXaxk27zwlk&feature=youtu.be&t=56m34s)
    unsigned int fastmod(const int input, const int ceil) {
      // apply the modulo operator only when needed
      return input >= ceil ? input % ceil : input;
    }
  }

  template<typename CharType>
  class WuManber {
  public:
    using StringType = std::basic_string<CharType>;

    WuManber(unsigned short HBITS = 4, size_t tableSize = 32768) :
      isInitialized_(false), m_(0), k_(0),
      HBITS_(HBITS), tableSize_(tableSize)
    {
      shiftTable_ = new size_t[tableSize_];
      hashPrefixTable_ = new std::vector<PatternHash>[tableSize_];
    }

    ~WuManber() {
      delete []shiftTable_;
      delete []hashPrefixTable_;
    }

    WuManber(const WuManber&) = delete;
    WuManber& operator =(const WuManber&) = delete;

    const std::vector<StringType>& patternList() const {
      return patternList_;
    }

    void preProcess(const std::vector<StringType> &patterns) {
      k_ = patterns.size();

      m_ = 0;
      for (const auto &pattern : patterns) {
        size_t patternLength = pattern.size();
        // throw an exception for length 1 and 2 patterns
        // those patterns can be handled by using something like lengthOneMap[ALPHABET_SIZE] and lengthTwoMap[ALPHABET_SIZE][ALPHABET_SIZE] for look up
        if (patternLength < B_) {
          throw std::logic_error("Pattern length must be at least 3");
        }
        m_ = (!m_) ? patternLength : std::min(patternLength, m_);
        patternList_.push_back(pattern);
      }

      // fill default value for SHIFT table
      for (int i = 0; i < tableSize_; ++i) {
        shiftTable_[i] = m_ - B_ + 1;
      }

      // fill HASH/PREFIX and SHIFT tables
      for (size_t i = 0; i < k_; ++i) {
        for (size_t j = m_; j >= B_; --j) {
          unsigned int hashValue;
          hashValue = patternList_[i][j - 1];
          hashValue <<= HBITS_;
          hashValue += patternList_[i][j - 1 - 1];
          hashValue <<= HBITS_;
          hashValue += patternList_[i][j - 2 - 1];
          hashValue = fastmod(hashValue, tableSize_);

          size_t shiftLength = m_ - j;
          shiftTable_[hashValue] = std::min(shiftTable_[hashValue], shiftLength);
          if (!shiftLength) {
            PatternHash patternHashToAdd;
            patternHashToAdd.idx = i;

            // calculate this prefixHash to help us skip some patterns if there are collisions in hashPrefixTable_
            patternHashToAdd.prefixHash = patternList_[i][0];
            patternHashToAdd.prefixHash <<= HBITS_;
            patternHashToAdd.prefixHash += patternList_[i][1];
            hashPrefixTable_[hashValue].push_back(patternHashToAdd);
          }
        }
      }

      isInitialized_ = true;
    }

    // onMatch takes 3 arguments: the matched pattern, the pattern's index in the pattern list, the start index of the match in text
    void scan(StringType text, std::function<void(const StringType&, size_t, size_t)> onMatch) {
      size_t textLength = text.size();
      if (!isInitialized_ || m_ > textLength) {
        return;
      }

      // index in text
      size_t idx = m_ - 1;
      while (idx < textLength) {
        // hash value for HASH table
        unsigned int hashValue;
        hashValue = text[idx];
        hashValue <<= HBITS_;
        hashValue += text[idx - 1];
        hashValue <<= HBITS_;
        hashValue += text[idx - 2];
        hashValue = fastmod(hashValue, tableSize_);

        size_t shiftLength = shiftTable_[hashValue];
        if (shiftLength == 0) {
          // found a potential match, check values in HASH/PREDIX and will shift by 1 character
          shiftLength = 1;

          // hash value to match pattern
          unsigned int prefixHash;
          prefixHash = text[idx - m_ + 1];
          prefixHash <<= HBITS_;
          prefixHash += text[idx - m_ + 2];
          for (const auto &potentialMatch : hashPrefixTable_[hashValue]) {
            if (prefixHash == potentialMatch.prefixHash) {
              bool isMatched = false;
              const StringType &pattern = patternList_[potentialMatch.idx];
              size_t idxInPattern = 0;
              size_t idxInText = idx - m_ + 1;
              size_t patternLength = pattern.size();

              // prefix hash matched so we try to match character by character
              while(idxInPattern < patternLength && idxInText < textLength && pattern[idxInPattern++] == text[idxInText++]);

              // end of pattern reached => match found
              if (idxInPattern == patternLength) {
                onMatch(pattern, potentialMatch.idx, idx - m_ + 1);
              }
            }
          }
        }
        idx += shiftLength;
      }

    }

  private:
    // block size
    // the paper says in practice, we use either B = 2 or B = 3
    // we'll use 3
    const size_t B_ = 3;

    // min pattern size
    size_t m_;

    // number of patterns
    size_t k_;

    // number of bits to shift when hashing
    // the paper says it use 5
    unsigned short HBITS_;

    // size of HASH and SHIFT tables
    size_t tableSize_;

    // SHIFT table
    size_t* shiftTable_;

    // store index in pattern list + prefix hash value for each pattern
    struct PatternHash
    {
        unsigned int prefixHash;
        size_t idx;
    };

    // HASH + PREFIX table
    std::vector<PatternHash>* hashPrefixTable_;

    // pattern list
    std::vector<StringType> patternList_;

    bool isInitialized_;
  };

} // namespace wu_manber

#endif // WU_MANBER_HPP

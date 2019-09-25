#ifndef WU_MANBER_HPP
#define WU_MANBER_HPP

#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <cmath>

namespace wu_manber {

  namespace { // anonymous namespace, things in here are "private" to wu_manber namespace

    // fast mod (ref: https://www.youtube.com/watch?v=nXaxk27zwlk&feature=youtu.be&t=56m34s)
    unsigned int fastmod(const int input, const int ceil) {
      // apply the modulo operator only when needed
      return input >= ceil ? input % ceil : input;
    }
  }

  template<typename CharType, typename CharTraits = std::char_traits<CharType>>
  class WuManber {
  public:
    using StringType = std::basic_string<CharType, CharTraits>;

    WuManber(unsigned short HBITS = 4, size_t tableSize = 32768) :
      isInitialized_(false), m_(0), k_(0),
      HBITS_(HBITS), tableSize_(tableSize)
    {
      shiftTable_ = new size_t[tableSize_];
      hashPrefixTable_ = new std::vector<PatternHash>[tableSize_];

      alphabetSize_ = pow(2, 8 * sizeof(CharType));
      isShortPatternExist_ = false;
    }

    ~WuManber() {
      delete []shiftTable_;
      delete []hashPrefixTable_;

      if (isShortPatternExist_) {
        delete []lengthOnePatternLookup_;
        for (int i = 0; i < alphabetSize_; ++i) {
          delete []lengthTwoPatternLookup_[i];
        }
        delete []lengthTwoPatternLookup_;
      }
    }

    WuManber(const WuManber&) = delete;
    WuManber& operator =(const WuManber&) = delete;

    const std::vector<StringType>& patternList() const {
      return patternList_;
    }

    void preProcess(const std::vector<StringType> &patterns) {
      m_ = 0;
      for (const auto &pattern : patterns) {
        size_t patternLength = pattern.size();
        if (patternLength < B_) {
          if (patternLength == 1) {
            lengthOnePatterns_.push_back(pattern);
          } else if (patternLength == 2) {
            lengthTwoPatterns_.push_back(pattern);
          }
          continue;
        }
        m_ = (!m_) ? patternLength : std::min(patternLength, m_);
        patternList_.push_back(pattern);
      }
      k_ = patternList_.size();

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

      isShortPatternExist_ = (lengthOnePatterns_.size() > 0) || (lengthTwoPatterns_.size() > 0);
      if (isShortPatternExist_) {
        lengthOnePatternLookup_ = new int[alphabetSize_];
        lengthTwoPatternLookup_ = new int*[alphabetSize_];
        for (int i = 0; i < alphabetSize_; ++i) {
          lengthOnePatternLookup_[i] = -1;
          lengthTwoPatternLookup_[i] = new int[alphabetSize_];
          for (int j = 0; j < alphabetSize_; ++j) {
            lengthTwoPatternLookup_[i][j] = -1;
          }
        }

        for (int i = 0; i < lengthOnePatterns_.size(); ++i) {
          lengthOnePatternLookup_[(size_t)lengthOnePatterns_[i][0]] = i;
        }

        for (int i = 0; i < lengthTwoPatterns_.size(); ++i) {
          lengthTwoPatternLookup_[(size_t)lengthTwoPatterns_[i][0]][(size_t)lengthTwoPatterns_[i][1]] = i;
        }
      }

      isInitialized_ = true;
    }

    // onMatch takes 3 arguments: the matched pattern, the pattern's index in the pattern list, the start index of the match in text
    void scan(const StringType &text, std::function<void(const StringType&, size_t, size_t)> onMatch) {
      size_t textLength = text.size();
      if (!isInitialized_ || textLength == 0) {
        return;
      }

      if (isShortPatternExist_) {
        int firstCharacterMatchIndex = lengthOnePatternLookup_[(size_t)text[0]];
        if (firstCharacterMatchIndex > -1) {
          onMatch(lengthOnePatterns_[firstCharacterMatchIndex], firstCharacterMatchIndex, 0);
        }
        const int PRE_WU_MANBER_LIMIT = std::min(m_ - 1, textLength);
        for (int idx = 1; idx < PRE_WU_MANBER_LIMIT; ++idx) {
          CharType preChar = text[idx - 1];
          CharType curChar = text[idx];
          checkShortPattern_(text, idx, onMatch);
        }
      }

      size_t idx = m_ - 1;
      while (idx < textLength) {
        if (isShortPatternExist_) {
          checkShortPattern_(text, idx, onMatch);
        }

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
        if (isShortPatternExist_) {
          ++idx;
        } else {
          idx += shiftLength;
        }
      }

    }

  private:
    // block size
    // the paper says in practice, we use either B = 2 or B = 3
    // we'll use 3
    const size_t B_ = 3;

    // min pattern size
    size_t m_;

    // number of patterns to be processed by Wu - Manber
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

    // handle length 1 and 2 patterns
    bool isShortPatternExist_;
    size_t alphabetSize_;
    int* lengthOnePatternLookup_;
    int** lengthTwoPatternLookup_;
    std::vector<StringType> lengthOnePatterns_;
    std::vector<StringType> lengthTwoPatterns_;

    bool isInitialized_;

    void checkShortPattern_(const StringType &text, size_t cur_idx, std::function<void(const StringType&, size_t, size_t)> onMatch) const {
      int l1MatchIndex = lengthOnePatternLookup_[(size_t)text[cur_idx]];
      if (l1MatchIndex > -1) {
        onMatch(lengthOnePatterns_[l1MatchIndex], l1MatchIndex, cur_idx);
      }

      int l2MatchIndex = lengthTwoPatternLookup_[(size_t)text[cur_idx - 1]][(size_t)text[cur_idx]];
      if (l2MatchIndex > -1) {
        onMatch(lengthTwoPatterns_[l2MatchIndex], l2MatchIndex, cur_idx - 1);
      }
    }
  };

} // namespace wu_manber

#endif // WU_MANBER_HPP

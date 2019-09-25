#ifndef WU_MANBER_HPP
#define WU_MANBER_HPP

#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <cstddef>
#include <climits>

namespace wu_manber {

  template<typename CharType, typename CharTraits = std::char_traits<CharType>>
  class WuManber {
  public:
    using StringType = std::basic_string<CharType, CharTraits>;

    WuManber(unsigned short int HBITS = 4, std::size_t tableSize = 32768) :
      isInitialized_(false), m_(0), k_(0),
      HBITS_(HBITS), tableSize_(tableSize),
      alphabetSize_(1 << (CHAR_BIT * sizeof(CharType))),
      isShortPatternExist_(false)
    {
      shiftTable_.resize(tableSize_);
      hashPrefixTable_.resize(tableSize_);
    }

    virtual ~WuManber() {}

    WuManber(const WuManber&) = delete;
    WuManber& operator =(const WuManber&) = delete;

    const std::vector<StringType>& patternList() const {
      return patternList_;
    }

    void preProcess(const std::vector<StringType> &patterns) {
      m_ = 0;
      for (const auto &pattern : patterns) {
        std::size_t patternLength = pattern.size();
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
      for (std::size_t i = 0; i < tableSize_; ++i) {
        shiftTable_[i] = m_ - B_ + 1;
      }

      // fill HASH/PREFIX and SHIFT tables
      for (std::size_t i = 0; i < k_; ++i) {
        for (std::size_t j = m_; j >= B_; --j) {
          unsigned int hashValue = getWuManberTableHashFromText_(patternList_[i], j - 1);

          std::size_t shiftLength = m_ - j;
          shiftTable_[hashValue] = std::min(shiftTable_[hashValue], shiftLength);
          if (!shiftLength) {
            PatternHash patternHashToAdd;
            patternHashToAdd.idx = i;

            // calculate this prefixHash to help us skip some patterns if there are collisions in hashPrefixTable_
            patternHashToAdd.prefixHash = getPrefixHashFromText_(patternList_[i], 0);
            hashPrefixTable_[hashValue].push_back(patternHashToAdd);
          }
        }
      }

      isShortPatternExist_ = (lengthOnePatterns_.size() > 0) || (lengthTwoPatterns_.size() > 0);
      if (isShortPatternExist_) {
        lengthOnePatternLookup_.resize(alphabetSize_);
        lengthTwoPatternLookup_.resize(alphabetSize_);
        for (std::size_t i = 0; i < alphabetSize_; ++i) {
          lengthOnePatternLookup_[i] = -1;
          lengthTwoPatternLookup_[i].resize(alphabetSize_);
          for (std::size_t j = 0; j < alphabetSize_; ++j) {
            lengthTwoPatternLookup_[i][j] = -1;
          }
        }

        for (std::size_t i = 0; i < lengthOnePatterns_.size(); ++i) {
          lengthOnePatternLookup_[static_cast<std::size_t>(lengthOnePatterns_[i][0])] = i;
        }

        for (std::size_t i = 0; i < lengthTwoPatterns_.size(); ++i) {
          lengthTwoPatternLookup_[static_cast<std::size_t>(lengthTwoPatterns_[i][0])][static_cast<std::size_t>(lengthTwoPatterns_[i][1])] = i;
        }
      }

      isInitialized_ = true;
    }

    // onMatch takes 3 arguments: the matched pattern, the pattern's index in the pattern list, the start index of the match in text
    void scan(const StringType &text, std::function<void(const StringType&, std::size_t, std::size_t)> onMatch) {
      std::size_t textLength = text.size();
      if (!isInitialized_ || textLength == 0) {
        return;
      }

      if (isShortPatternExist_) {
        int firstCharacterMatchIndex = lengthOnePatternLookup_[static_cast<std::size_t>(text[0])];
        if (firstCharacterMatchIndex > -1) {
          onMatch(lengthOnePatterns_[firstCharacterMatchIndex], firstCharacterMatchIndex, 0);
        }
        const int PRE_WU_MANBER_LIMIT = std::min(m_ - 1, textLength);
        for (std::size_t idx = 1; idx < PRE_WU_MANBER_LIMIT; ++idx) {
          CharType preChar = text[idx - 1];
          CharType curChar = text[idx];
          checkShortPattern_(text, idx, onMatch);
        }
      }

      std::size_t idx = m_ - 1;
      while (idx < textLength) {
        if (isShortPatternExist_) {
          checkShortPattern_(text, idx, onMatch);
        }

        // hash value for HASH table
        unsigned int hashValue = getWuManberTableHashFromText_(text, idx);

        std::size_t shiftLength = shiftTable_[hashValue];
        if (shiftLength == 0) {
          // found a potential match, check values in HASH/PREDIX and will shift by 1 character
          shiftLength = 1;

          // hash value to match pattern
          unsigned int prefixHash = getPrefixHashFromText_(text, idx - m_ + 1);
          for (const auto &potentialMatch : hashPrefixTable_[hashValue]) {
            if (prefixHash == potentialMatch.prefixHash) {
              bool isMatched = false;
              const StringType &pattern = patternList_[potentialMatch.idx];
              std::size_t idxInPattern = 0;
              std::size_t idxInText = idx - m_ + 1;
              std::size_t patternLength = pattern.size();

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
    static constexpr std::size_t B_ = 3;

    // min pattern size
    std::size_t m_;

    // number of patterns to be processed by Wu - Manber
    std::size_t k_;

    // number of bits to shift when hashing
    // the paper says it use 5
    unsigned short int HBITS_;

    // size of HASH and SHIFT tables
    std::size_t tableSize_;

    // SHIFT table
    std::vector<std::size_t> shiftTable_;

    // store index in pattern list + prefix hash value for each pattern
    struct PatternHash
    {
        unsigned int prefixHash;
        std::size_t idx;
    };

    // HASH + PREFIX table
    std::vector<std::vector<PatternHash>> hashPrefixTable_;

    // pattern list
    std::vector<StringType> patternList_;

    // handle length 1 and 2 patterns
    bool isShortPatternExist_;
    std::size_t alphabetSize_;
    std::vector<std::size_t> lengthOnePatternLookup_;
    std::vector<std::vector<std::size_t>> lengthTwoPatternLookup_;
    std::vector<StringType> lengthOnePatterns_;
    std::vector<StringType> lengthTwoPatterns_;

    bool isInitialized_;

    void checkShortPattern_(const StringType &text, std::size_t cur_idx, std::function<void(const StringType&, std::size_t, std::size_t)> onMatch) const {
      int l1MatchIndex = lengthOnePatternLookup_[static_cast<std::size_t>(text[cur_idx])];
      if (l1MatchIndex > -1) {
        onMatch(lengthOnePatterns_[l1MatchIndex], l1MatchIndex, cur_idx);
      }

      int l2MatchIndex = lengthTwoPatternLookup_[static_cast<std::size_t>(text[cur_idx - 1])][static_cast<std::size_t>(text[cur_idx])];
      if (l2MatchIndex > -1) {
        onMatch(lengthTwoPatterns_[l2MatchIndex], l2MatchIndex, cur_idx - 1);
      }
    }

    unsigned int getWuManberTableHashFromText_(const StringType &text, std::size_t cur_idx) {
      unsigned int hashValue;
      hashValue = text[cur_idx];
      hashValue <<= HBITS_;
      hashValue += text[cur_idx - 1];
      hashValue <<= HBITS_;
      hashValue += text[cur_idx - 2];
      hashValue %= tableSize_; // may try this for faster mod: https://www.youtube.com/watch?v=nXaxk27zwlk&feature=youtu.be&t=56m34s
      return hashValue;
    }

    unsigned int getPrefixHashFromText_(const StringType &text, std::size_t cur_idx) {
      unsigned int prefixHash;
      prefixHash = text[cur_idx];
      prefixHash <<= HBITS_;
      prefixHash += text[cur_idx + 1];
      return prefixHash;
    }
  };

} // namespace wu_manber

#endif // WU_MANBER_HPP

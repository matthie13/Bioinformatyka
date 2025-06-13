#include "levenshtein.h"

int levenshtein_score(const std::string& s, const std::string& t) {
    if (s == t) return 0;

    const size_t n = s.length();
    const size_t m = t.length();

    if (n == 0) return static_cast<int>(m);
    if (m == 0) return static_cast<int>(n);

    std::vector<int> p(n + 1);
    std::vector<int> d(n + 1);

    for (size_t i = 0; i <= n; i++) {
        p[i] = static_cast<int>(i);
    }

    for (size_t j = 0; j <= m; j++) {
        char tJ = t[j - 1];
        d[0] = static_cast<int>(j);

        for (size_t i = 1; i <= n; i++) {
            int cost = s[i - 1] == tJ ? 0 : 1;
            d[i] = std::min(std::min(d[i - 1] + 1, p[i] + 1), p[i - 1] + cost);
        }
        std::swap(p, d);
    }
    return p[n];
}

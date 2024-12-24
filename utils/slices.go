package utils

func MaxIndexFunc(N int, gt func(int, int) bool) int {
	res := 0
	for i := 1; i < N; i++ {
		if gt(i, res) {
			res = i
		}
	}
	return res
}

func PartialSums(s ...int) []int {
	if len(s) == 0 {
		return nil
	}
	sums := make([]int, len(s))
	sums[0] = s[0]
	for i := 1; i < len(s); i++ {
		sums[i] = sums[i-1] + s[i]
	}
	return sums
}

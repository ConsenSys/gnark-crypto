//go:build !forcegen

package git

import (
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
)

var gitStatusDiff string
var once sync.Once

func HasChanges(dir string) bool {
	dir = filepath.Clean(dir)
	once.Do(func() {
		cmd := exec.Command("git", "status", "--porcelain")
		output, _ := cmd.Output()
		gitStatusDiff = string(output)
	})

	// if status contains bavard or addchain, we return true by default
	if strings.Contains(gitStatusDiff, "go.mod") || strings.Contains(gitStatusDiff, "go.sum") {
		return true
	}

	return strings.Contains(gitStatusDiff, dir)
}

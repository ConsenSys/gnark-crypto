package pairing

import (
	"path/filepath"
	"strings"

	"github.com/consensys/bavard"
	"github.com/consensys/gnark-crypto/internal/generator/config"
	"github.com/consensys/gnark-crypto/internal/generator/git"
)

func Generate(conf config.Curve, baseDir string, bgen *bavard.BatchGenerator) error {
	if !(git.HasChanges("./pairing/template")) {
		return nil
	}
	packageName := strings.ReplaceAll(conf.Name, "-", "")
	return bgen.Generate(conf, packageName, "./pairing/template", bavard.Entry{
		File: filepath.Join(baseDir, "pairing_test.go"), Templates: []string{"tests/pairing.go.tmpl"},
	})

}

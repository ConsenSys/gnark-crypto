package pedersen

import (
	"path/filepath"

	"github.com/consensys/bavard"
	"github.com/consensys/gnark-crypto/internal/generator/config"
)

func Generate(conf config.Curve, baseDir string, bgen *bavard.BatchGenerator) error {
	// pedersen commitment scheme
	conf.Package = "pedersen"
	entries := []bavard.Entry{
		{File: filepath.Join(baseDir, "doc.go"), Templates: []string{"doc.go.tmpl"}},
		{File: filepath.Join(baseDir, "pedersen.go"), Templates: []string{"pedersen.go.tmpl"}},
		{File: filepath.Join(baseDir, "pedersen_test.go"), Templates: []string{"pedersen.test.go.tmpl"}},
		{File: filepath.Join(baseDir, "example_test.go"), Templates: []string{"example_test.go.tmpl"}},
	}
	return bgen.Generate(conf, conf.Package, "./pedersen/template/", entries...)

}

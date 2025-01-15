package fflonk

import (
	"path/filepath"

	"github.com/consensys/bavard"
	"github.com/consensys/gnark-crypto/internal/generator/config"
	"github.com/consensys/gnark-crypto/internal/generator/git"
)

func Generate(conf config.Curve, baseDir string, bgen *bavard.BatchGenerator) error {
	if !(git.HasChanges("./fflonk/template")) {
		return nil
	}
	// kzg commitment scheme
	conf.Package = "fflonk"
	entries := []bavard.Entry{
		{File: filepath.Join(baseDir, "doc.go"), Templates: []string{"doc.go.tmpl"}},
		{File: filepath.Join(baseDir, "fflonk.go"), Templates: []string{"fflonk.go.tmpl"}},
		{File: filepath.Join(baseDir, "fflonk_test.go"), Templates: []string{"fflonk.test.go.tmpl"}},
		{File: filepath.Join(baseDir, "marshal.go"), Templates: []string{"marshal.go.tmpl"}},
		{File: filepath.Join(baseDir, "example_test.go"), Templates: []string{"example_test.go.tmpl"}},
	}
	return bgen.Generate(conf, conf.Package, "./fflonk/template/", entries...)

}

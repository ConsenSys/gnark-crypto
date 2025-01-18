package mpcsetup

import (
	"github.com/consensys/bavard"
	"github.com/consensys/gnark-crypto/internal/generator/config"
	"path/filepath"
)

func Generate(conf config.Curve, baseDir string, bgen *bavard.BatchGenerator) error {
	entries := []bavard.Entry{
		{File: filepath.Join(baseDir, "doc.go"), Templates: []string{"doc.go.tmpl"}},
		{File: filepath.Join(baseDir, "mpcsetup.go"), Templates: []string{"mpcsetup.go.tmpl"}},
		{File: filepath.Join(baseDir, "mpcsetup_test.go"), Templates: []string{"tests/mpcsetup.go.tmpl"}},
	}
	return bgen.Generate(conf, "mpcsetup", "./mpcsetup/template", entries...)
}

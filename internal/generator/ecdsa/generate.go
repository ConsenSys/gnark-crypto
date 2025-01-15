package ecdsa

import (
	"path/filepath"

	"github.com/consensys/bavard"
	"github.com/consensys/gnark-crypto/internal/generator/config"
	"github.com/consensys/gnark-crypto/internal/generator/git"
)

func Generate(conf config.Curve, baseDir string, bgen *bavard.BatchGenerator) error {
	if !(git.HasChanges("./ecdsa/template")) {
		return nil
	}
	// ecdsa
	conf.Package = "ecdsa"
	baseDir = filepath.Join(baseDir, conf.Package)

	entries := []bavard.Entry{
		{File: filepath.Join(baseDir, "doc.go"), Templates: []string{"doc.go.tmpl"}},
		{File: filepath.Join(baseDir, "ecdsa.go"), Templates: []string{"ecdsa.go.tmpl"}},
		{File: filepath.Join(baseDir, "ecdsa_test.go"), Templates: []string{"ecdsa.test.go.tmpl"}},
		{File: filepath.Join(baseDir, "marshal.go"), Templates: []string{"marshal.go.tmpl"}},
		{File: filepath.Join(baseDir, "marshal_test.go"), Templates: []string{"marshal.test.go.tmpl"}},
	}
	return bgen.Generate(conf, conf.Package, "./ecdsa/template", entries...)

}

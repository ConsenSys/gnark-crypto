package hash_to_field

import (
	"path/filepath"

	"github.com/consensys/bavard"
	"github.com/consensys/gnark-crypto/internal/generator/config"
)

func Generate(conf config.FieldDependency, baseDir string, bgen *bavard.BatchGenerator) error {
	entries := []bavard.Entry{
		{File: filepath.Join(baseDir, "doc.go"), Templates: []string{"doc.go.tmpl"}},
		{File: filepath.Join(baseDir, "hash_to_field.go"), Templates: []string{"hash_to_field.go.tmpl"}},
		{File: filepath.Join(baseDir, "hash_to_field_test.go"), Templates: []string{"hash_to_field_test.go.tmpl"}},
	}

	return bgen.Generate(conf, "hash_to_field", "./hash_to_field/template", entries...)
}

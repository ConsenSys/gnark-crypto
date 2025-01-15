package main

import (
	"encoding/json"
	"fmt"
	"os"

	"github.com/consensys/gnark-crypto/ecc/bn254/fr"
	"github.com/consensys/gnark-crypto/ecc/bn254/fr/mimc"
	"github.com/consensys/gnark-crypto/internal/generator/git"
)

type numericalMiMCTestCase struct {
	In  []string `json:"in"`
	Out string   `json:"out"`
}

func assertNoError(err error) {
	if err != nil {
		fmt.Println(err)
		os.Exit(-1)
	}
}

//go:generate go run main.go
func main() {
	if !git.HasChanges("mimc/vectors.json") {
		fmt.Println("no changes in mimc/vectors.json, skipping generation")
		return
	}
	fmt.Println("generating test vectors for MiMC...")
	var tests []numericalMiMCTestCase

	bytes, err := os.ReadFile("./vectors.json")
	assertNoError(err)
	assertNoError(json.Unmarshal(bytes, &tests))

	hsh := mimc.NewMiMC()

	for i := range tests {

		hsh.Reset()
		var x fr.Element
		for j := range tests[i].In {
			_, err = x.SetString(tests[i].In[j])
			assertNoError(err)

			b := x.Bytes()
			_, err = hsh.Write(b[:])
			assertNoError(err)
		}

		bytes = hsh.Sum(nil)

		x.SetBytes(bytes)
		tests[i].Out = "0x" + x.Text(16)
	}

	bytes, err = json.MarshalIndent(tests, "", "\t")
	assertNoError(err)
	err = os.WriteFile("./vectors.json", bytes, 0)
	assertNoError(err)
}

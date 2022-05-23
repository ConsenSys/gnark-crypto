package element

const OpsNoAsm = `
// /!\ WARNING /!\
// this code has not been audited and is provided as-is. In particular, 
// there is no security guarantees such as constant time implementation 
// or side-channel attack resistance
// /!\ WARNING /!\

{{ $mulConsts := list 3 5 13 }}
{{- range $i := $mulConsts }}

// MulBy{{$i}} x *= {{$i}}
func MulBy{{$i}}(x *{{$.ElementName}}) {
	{{- if eq 1 $.NbWords}}
	var y {{$.ElementName}}
	y.SetUint64({{$i}})
	x.Mul(x, &y)
	{{- else}}
	mulByConstant(x, {{$i}})
	{{- end}}
}

{{- end}}


// Butterfly sets 
// a = a + b
// b = a - b 
func Butterfly(a, b *{{.ElementName}}) {
	_butterflyGeneric(a, b)
}

{{- if ne .NbWords 1}}
func mul(z, x, y *{{.ElementName}}) {
	_mulGeneric(z, x, y)
}
{{- end}}


// FromMont converts z in place (i.e. mutates) from Montgomery to regular representation
// sets and returns z = z * 1
func fromMont(z *{{.ElementName}} ) {
	_fromMontGeneric(z)
}


func neg(z,  x *{{.ElementName}}) {
	_negGeneric(z,x)
}

func reduce(z *{{.ElementName}})  {
	_reduceGeneric(z)
}
`

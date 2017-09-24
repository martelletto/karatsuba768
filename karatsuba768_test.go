// Copyright (c) 2017 Pedro Martelletto. All rights reserved.
// Use of this source code is governed by a BSD-style license
// that can be found in the LICENSE file.

package karatsuba768

import (
	"bufio"
	"compress/gzip"
	"errors"
	"fmt"
	"io"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"testing"
)

func TestFreeze(t *testing.T) {
	for i := int32(-165191049); i < 165191050; i++ {
		x := Freeze(i)
		y := ((i % 9829) + 9829) % 9829
		if x != y {
			t.Fatalf("x=%d != y=%d for i=%d", x, y, i)
		}
	}
}

func loadPoly(buf *bufio.Reader, p []int32, size int) error {
	body, chunk, err := buf.ReadLine()
	if chunk {
		return errors.New("input buffer too small")
	}
	if err != nil {
		return err
	}
	line := string(body)
	parts := strings.Split(line, ",")
	if len(parts) > size {
		return errors.New("too many parts")
	}
	for i := range parts {
		n, err := strconv.ParseInt(strings.TrimSpace(parts[i]), 10, 32)
		if err != nil {
			return err
		}
		p[i] = int32(n)
	}
	return nil
}

func cmpPoly(t *testing.T, c, d *[1536]int32) error {
	for i := 0; i < 1536; i++ {
		if c[i] != d[i] {
			return fmt.Errorf("c=%d, d=%d for i=%d", c[i], d[i], i)
		}
	}
	return nil
}


func textbookMul(h *[1536]int32, f, g *[768]int32) {
	for i := 0; i < 768; i++ {
		for j := 0; j < 768; j++ {
			h[i+j] += (f[i] * g[j]) % 9829
		}
	}
	for i := 0; i < 1536; i++ {
		h[i] %= 9829
	}
}

func TestSage64(t *testing.T) {
	// open file
	f, err := os.Open("sage64.gz")
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()

	// create input buf
	in, err := gzip.NewReader(f)
	if err != nil {
		t.Fatal(err)
	}
	defer in.Close()
	buf := bufio.NewReaderSize(in, 1 << 14)

	for ln := 0;; ln += 3 {
		// load a
		t.Logf("processing line %d", ln + 1)
		a := new([768]int32)
		err := loadPoly(buf, a[:], 768)
		if err != nil {
			if err == io.EOF {
				break
			}
			t.Fatalf("couldn't read poly: %v", err)
		}

		// load b
		t.Logf("processing line %d", ln + 2)
		b := new([768]int32)
		err = loadPoly(buf, b[:], 768)
		if err != nil {
			t.Fatalf("couldn't read poly: %v", err)
		}

		// load c
		t.Logf("processing line %d", ln + 3)
		c := new([1536]int32)
		err = loadPoly(buf, c[:], 1536)
		if err != nil {
			t.Fatalf("couldn't read poly: %v", err)
		}

		// verify that c == d
		d := new([1536]int32)
		Mul(d, a, b)
		err = cmpPoly(t, c, d)
		if err != nil {
			t.Fatalf("c != d: %v", err)
		}
	}
}

func TestRandom64(t *testing.T) {
	for i := 0; i < 64; i++ {
		a := new([768]int32)
		b := new([768]int32)
		for j := 0; j < 768; j++ {
			a[j] = int32(rand.Intn(9829))
			b[j] = int32(rand.Intn(9829))
		}
		c := new([1536]int32)
		d := new([1536]int32)
		textbookMul(c, a, b)
		Mul(d, a, b)
		err := cmpPoly(t, c, d)
		if err != nil {
			fmt.Fprintf(os.Stderr, "a=%v\n", a)
			fmt.Fprintf(os.Stderr, "b=%v\n", b)
			fmt.Fprintf(os.Stderr, "c=%v\n", c)
			fmt.Fprintf(os.Stderr, "d=%v\n", d)
			t.Fatalf("c != d: %v", err)
		}
	}
}

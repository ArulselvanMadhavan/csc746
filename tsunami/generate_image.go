package main

import (
	"image"
	"image/color"
	"image/draw"
	"image/png"
	"log"
	"os"
	"path/filepath"
	"strings"
	"sync"

	"github.com/lucasb-eyer/go-colorful"
)

const GRAPHICS_PATH = "/graphics_output/*.data"

const WIDTH = 800
const HEIGHT = 800

func processFile(wg *sync.WaitGroup, fpath string) {
	defer wg.Done()
	fname := filepath.Base(fpath)
	fname = strings.Split(fname, ".")[0]
	upLeft := image.Point{0, 0}
	lowRight := image.Point{WIDTH, HEIGHT}

	img := image.NewRGBA(image.Rectangle{upLeft, lowRight})
	fdir := filepath.Dir(fpath)
	iFile, err := os.Create(fdir + "/" + fname + ".png")
	if err != nil {
		log.Println("Unable to create file", err)
	}
	s_color := colorful.Hsl(216.0, 1.0, 0.5)
	red_rect := image.Rect(0, 0, 120, 160) //  geometry of 2nd rectangle

	r, g, b, a := s_color.RGBA()
	myred := color.RGBA{uint8(r), uint8(g), uint8(b), uint8(a)}
	draw.Draw(img, red_rect, &image.Uniform{myred}, image.ZP, draw.Src)
	// img.Rect
	png.Encode(iFile, img)
	log.Println(fpath)
}

func main() {
	var wg sync.WaitGroup
	path, _ := os.Getwd()
	if strings.Contains(path, "build") {
		path += GRAPHICS_PATH
	} else {
		path += "/build" + GRAPHICS_PATH
	}
	log.Println("Path:", path)
	matches, err := filepath.Glob(path)
	if err != nil {
		log.Fatalf("error reading .data files. %v", err)
	}

	for _, fi := range matches {
		wg.Add(1)
		go processFile(&wg, fi)
		// log.Println(strings.Split(fi, " "))
	}
	wg.Wait()
	log.Println("All done")
}

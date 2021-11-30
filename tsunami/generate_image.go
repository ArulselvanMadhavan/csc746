package main

import (
	"bufio"
	"bytes"
	"fmt"
	"image"
	"image/color"
	"image/draw"
	"image/png"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"sync"

	"github.com/lucasb-eyer/go-colorful"
	"golang.org/x/image/font"
	"golang.org/x/image/font/inconsolata"
	"golang.org/x/image/math/fixed"
	"gonum.org/v1/hdf5"
)

const GRAPHICS_PATH = "/graphics_output/*.data"

const WIDTH = 800
const HEIGHT = 800

func addLabel(img *image.RGBA, x, y int, label string) {
	col := color.RGBA{255, 255, 255, 255}
	point := fixed.Point26_6{fixed.Int26_6(x * 128), fixed.Int26_6(y * 128)}
	fce := inconsolata.Bold8x16
	fce.Height = fce.Height * 2
	fce.Width = fce.Width * 2
	d := &font.Drawer{
		Dst:  img,
		Src:  image.NewUniform(col),
		Face: fce,
		Dot:  point,
	}
	d.DrawString(label)
}

func buildPixel(chunks [][]byte, points []int) (image.Rectangle, color.RGBA) {
	for i, c := range chunks {
		points[i], _ = strconv.Atoi(string(c))
	}
	x1 := points[0]
	x2 := points[0] + points[2]
	y1 := points[1]
	y2 := points[1] + points[3]
	p_color := points[4]
	// n_color := float64(p_color) / 255.0
	// col := colorful.Color{n_color, n_color, n_color}
	// colorful.FastLinearRgb(r float64, g float64, b float64)
	col := colorful.Hsl(float64(p_color), 1.0, 0.5)
	r, g, b, a := col.RGBA()
	rgba_col := color.RGBA{uint8(r), uint8(g), uint8(b), uint8(a)}
	pixel := image.Rect(x1, y1, x2, y2)
	return pixel, rgba_col
}

func buildLabel(chunks [][]byte) (int, float64) {
	iteration, _ := strconv.Atoi(string(chunks[0]))
	time, err := strconv.ParseFloat(string(chunks[1]), 64)
	if err != nil {
		log.Fatalf("Float conversion failed. %s", err)
	}
	return iteration, time
}

func processFile(wg *sync.WaitGroup, fpath string) {
	defer wg.Done()
	points := make([]int, 5)
	fname := getFileName(fpath)
	// Create image
	upLeft := image.Point{0, 0}
	lowRight := image.Point{WIDTH, HEIGHT}

	img := image.NewRGBA(image.Rectangle{upLeft, lowRight})

	// Read file
	f, err := os.Open(fpath)
	defer f.Close()
	if err != nil {
		log.Fatalf("File open failed. Unable to read .data file. %s", err)
	}
	reader := bufio.NewReader(f)
	var iteration int
	var simTime float64
	for {
		line, err := reader.ReadBytes('\n')
		if len(line) == 0 {
			break
		}
		if err != nil {
			log.Fatalf("Error occured when reading file:%s\nError:%s", fpath, err)
		}
		line = bytes.TrimRight(line, "\n")
		chunks := bytes.Split(line, []byte(","))
		if len(chunks) == 2 {
			iteration, simTime = buildLabel(chunks)
		} else if len(chunks) == 5 {
			pixel, rgba_col := buildPixel(chunks, points)
			draw.Draw(img, pixel, &image.Uniform{rgba_col}, image.ZP, draw.Src)
		} else {
			log.Fatalf("Line contains unexpected number of chunks", len(chunks), line)
		}
	}

	// Stick labels
	addLabel(img, 0, 13, "Iteration: "+strconv.Itoa(iteration))
	addLabel(img, 0, 26, "time: "+fmt.Sprint(simTime))

	fdir := filepath.Dir(fpath)
	iFile, err := os.Create(fdir + "/" + fname + ".png")
	if err != nil {
		log.Println("Unable to create file", err)
	}
	png.Encode(iFile, img)
}

func getFileName(fpath string) string {
	return strings.Split(filepath.Base(fpath), ".")[0]
}

func createImageList(matches []string) {
	f, err := os.Create(filepath.Dir(matches[0]) + "/imagelist.txt")
	defer f.Close()
	if err != nil {
		log.Fatal("Unable to create imagelist file")
	}
	writer := bufio.NewWriter(f)
	for _, fpath := range matches {
		if writer.Size() > 0 {
			writer.WriteString("\n")
		}
		fname := getFileName(fpath) + ".png"
		writer.WriteString(fname)
	}
	writer.Flush()
}

func main2() {
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
	}
	wg.Wait()
	// Create image list
	createImageList(matches)
	log.Println("All done")
}

func main() {
	fname := "example0.hdf5"
	f, err := hdf5.OpenFile(fname, hdf5.F_ACC_RDONLY)
	if err != nil {
		panic(err)
	}
	dset, err := f.OpenDataset("data array")
	if err != nil {
		panic(err)
	}
	dspace := dset.Space()
	dims, _, _ := dspace.SimpleExtentDims()
	totalDims := 1
	for _, v := range dims {
		totalDims = totalDims * int(v)
	}
	s2 := make([]float64, totalDims)
	err = dset.Read(&s2)
	log.Println("Size:", totalDims)
	log.Println("All done", s2)
}

package main

import (
	"fmt"
	"flag"
	"os"
	"bufio"
	"strings"
	"strconv"
	"launchpad.net/mgo"
)

type Location struct {
	FileName    string
	Start       int64
	End         int64
	Strand      string
	Prob        float64
	NormScore   float64
	AggScore    float64
	Matrix      string
	Subsequence string
	Pval        float64
}

func LocationFromLine(line []byte) *Location {
	vals := strings.Split(string(line), "\t")
	return &Location{FileName: vals[0],
		Start:       parseInt64(vals[1]),
		End:         parseInt64(vals[2]),
		Strand:      vals[3],
		Prob:        parseFloat64(vals[4]),
		NormScore:   parseFloat64(vals[5]),
		AggScore:    parseFloat64(vals[6]),
		Matrix:      vals[7],
		Subsequence: vals[8],
		Pval:        parseFloat64(vals[9])}
}

func parseInt64(s string) int64 {
	i, err := strconv.Atoi64(s)
	if err != nil {
		panic(err)
	}
	return i
}

func parseFloat64(s string) float64 {
	i, err := strconv.Atof64(s)
	if err != nil {
		panic(err)
	}
	return i
}

func main() {
	//parse args
	var file string
	var mongo string
	var dbname string
	var collection string

	flag.StringVar(&file, "file", "", "The file to load.")
	flag.StringVar(&mongo, "mongo", "", "The mongo server to connect to.")
	flag.StringVar(&dbname, "db", "molotov", "The bb name to use.")
	flag.StringVar(&collection, "c", "locations", "The name of the collection to load data into.")
	flag.Parse()

	fmt.Printf("Connection to mongo.\n")
	//open connection to mongo
	session, err := mgo.Mongo(mongo)
	if err != nil {
		panic(err)
	}
	defer session.Close()

	// Optional. Switch the session to a monotonic behavior.
	//session.SetMode(mgo.Monotonic, true)

	fmt.Printf("Opening collection.\n")
	//select db and collection to use
	c := session.DB(dbname).C(collection)

	fmt.Printf("Opening file.\n")
	//open file
	fo, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer fo.Close()

	reader := bufio.NewReader(fo)
	for {
		line, ispre, err := reader.ReadLine()
		if ispre {
			panic("Line was too long.")
		}
		if err == os.EOF {
			break
		}
		if err != nil {
			panic(err)
		}

		err = c.Insert(LocationFromLine(line))
		if err != nil {
			panic(err)
		}

	}

}

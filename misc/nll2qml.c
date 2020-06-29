/*
 * Copyright (C) 1999-2010 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser Public License for more details.

 * You should have received a copy of the GNU Lesser Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */


/**   nll2qml.c

	Program to convert NLL Hypocenter-Phase format to QuakeML

 */

/**-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/**
	history:	(see also http://alomax.net/nlloc -> Updates)

	ver 01    22DEC2006  AJL  Original version developed from testWriter.c

 */


/** References */

/**
 * section: xmlWriter
 * synopsis: use various APIs for the xmlWriter
 * purpose: tests a number of APIs for the xmlWriter, especially
 *          the various methods to write to a filename, to a memory
 *          buffer, to a new document, or to a subtree. It shows how to
 *          do encoding string conversions too. The resulting
 *          documents are then serialized.
 * usage: testWriter
 * test: testWriter ; for i in 1 2 3 4 ; do diff writer.xml writer$$i.res ; done ; rm writer*.res
 * author: Alfred Mickautsch

 * copy: see Copyright for the status of this software.

The MIT License

Copyright (c) <year> <copyright holders>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Copyright © 2006 by the Open Source Initiative
The contents of this website are licensed under the Open Software License 2.1 or Academic Free License 2.1

OSI is a registered non-profit with 501(c)(3) status. Donating to OSI is one way to show your support.

 */


/** QuakeML example (from http://www.quakeml.ethz.ch/quakeml)

<?xml version="1.0" encoding="ISO-8859-1"?>
<location main="true" unique_id="LOC2004-03-18T11h16m55.651s" xmlns="http://quakeml.ethz.ch/ns/eventlist/location">
  <origin-time timezone="00:00" xmlns="http://quakeml.ethz.ch/ns/eventlist/location">
    <year>2004</year>
    <month>03</month>
    <day>17</day>
    <hour>16</hour>
    <minute>03</minute>
    <second>10.3</second>
  </origin-time>
  <latitude unit="degree">47.1</latitude>
  <longitude unit="degree">9.1</longitude>
  <depth unit="km">0</depth>
  <magnitude unit="ML">1.5</magnitude>
  <region>Kloental / Switzerland</region>
</location>

 */



#include <stdio.h>
#include <string.h>
#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>

#include "../src/GridLib.h"

// libxml, check
#if defined(LIBXML_WRITER_ENABLED) && defined(LIBXML_OUTPUT_ENABLED)

// libxml
//#define MY_ENCODING "ISO-8859-1"
#define MY_ENCODING "utf-8"
int initLibXml(void);
int cleanUpLibXml(void);
xmlChar *ConvertInput(const char *in, const char *encoding);

// QuakeML
xmlTextWriterPtr startQuakeMLDocument(const char *xmlWriterUri);
void endQuakeMLDocument(xmlTextWriterPtr writer);
void startQuakeMLEventElement(xmlTextWriterPtr writer, const char* uniqueId);

// nll
int nllHyp2Qml(const char *xmlWriterUri, const char *hyp_file_name);

// nll QuakeML
void writeQuakeMLLocationNLL(xmlTextWriterPtr writer,
			  const char* uniqueId, const char* timezone,
			  HypoDesc* phypo, ArrivalDesc* parrivals, int narrivals,
			  int iWriteArrivals, int iWriteMinimal,
			  GridDesc* pgrid, int n_proj);
void writeQuakeMLOriginTime(xmlTextWriterPtr writer, const char* timezone,
			    int year, int month, int day, int hour, int min, double sec);
void writeQuakeMLLatitude(xmlTextWriterPtr writer, double latitude, const char* units, double error);
void writeQuakeMLLongitude(xmlTextWriterPtr writer, double longitude, const char* units, double error);
void writeQuakeMLDepth(xmlTextWriterPtr writer, double depth, const char* units, double error);
void writeQuakeMLMagnitude(xmlTextWriterPtr writer, double magnitude, const char* units, double error);
void writeQuakeMLRegion(xmlTextWriterPtr writer, const char* region);



/** program to convert NLL hypo to QuakeML */

#define PNAME  "nll2qml"


int main(int argc, char** argv)
{

	int istat, narg;


	/* set program name */

	strcpy(prog_name, PNAME);


	/* check command line for correct usage */

	fprintf(stdout, "%s Arguments: ", prog_name);
	for (narg = 0; narg < argc; narg++)
		fprintf(stdout, "<%s> ", argv[narg]);
	fprintf(stdout, "\n");

	if (argc < 3) {
		nll_puterr("ERROR wrong number of command line arguments.");
		disp_usage(prog_name,
			   "<input_hyp_file> <output_qml_file_root_path>");
		exit(-1);
	}

	char hyp_file_name[FILENAME_MAX];
	char qml_file_name[FILENAME_MAX];
	strcpy(hyp_file_name, argv[1]);
	strcpy(qml_file_name, argv[2]);
	char* chr = strrchr(hyp_file_name, '/');
	if (chr == NULL)
		chr = hyp_file_name;
	else
		chr++;
	strcat(qml_file_name, chr);
	strcat(qml_file_name, ".xml");
	fprintf(stdout, "output_qml_file_name <%s>\n", qml_file_name);



	SetConstants();
	message_flag = 1;
	DispProgInfo();
	message_flag = 0;


	if ((istat =  nllHyp2Qml(qml_file_name, hyp_file_name)) < 0) {
		nll_puterr("ERROR doing nll2qml process.");
		exit(-1);
	}



	exit(0);

}



int nllHyp2Qml(const char *xmlWriterUri, const char *hyp_file_name)
{
	// this initialize the library
	initLibXml();

	// open new QuakeML document
	xmlTextWriterPtr writer = startQuakeMLDocument(xmlWriterUri);

	// open NLL hypocenter file
	HypoDesc Hypo;
	GridDesc locgrid;
	int istat = GetHypLoc(NULL, hyp_file_name, &Hypo, Arrival, &NumArrivals, 1, &locgrid, 0);
	if (istat < 0) {
		exit(-1);
	}

	// start new QuakeML event
	char uniqueId1[32];
	HypoDesc* phypo = &Hypo;
	sprintf(uniqueId1, "%4.4d%2.2d%2.2d.%2.2d%2.2d%5.2lf",
		phypo->year, phypo->month, phypo->day,
		phypo->hour, phypo->min, (double) phypo->sec);
	startQuakeMLEventElement(writer, uniqueId1);

	int iWriteArrivals = 0;
	int iWriteMinimal = 0;
	int n_proj = 0;

	char uniqueId2[32];
	strcpy(uniqueId2, uniqueId1);
	strcat(uniqueId2, "A");

	char* timezone = "00:00";

	// write QuakeML location
	writeQuakeMLLocationNLL(writer, uniqueId2, timezone, &Hypo, Arrival, NumArrivals,
			     iWriteArrivals, iWriteMinimal, &locgrid, n_proj);


	// end QuakeML document
	endQuakeMLDocument(writer);

	// Cleanup function for the XML library.
	cleanUpLibXml();

	/*
	* this is to debug memory for regression tests
	*/
	xmlMemoryDump();

	return 0;
}



/**
 * startQuakeMLDocument:
 * @uri: the output URI
 *
 */
xmlTextWriterPtr startQuakeMLDocument(const char *xmlWriterUri)
{
	int rc;
	xmlTextWriterPtr writer;

    // Create a new XmlWriter for xmlWriterUri, with no compression.
	writer = xmlNewTextWriterFilename(xmlWriterUri, 0);
	if (writer == NULL) {
		printf("nll2qml: Error creating the xml writer\n");
		return(writer);
	}

    // Set indenting on
	rc = xmlTextWriterSetIndent(writer, 1);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterSetIndent\n");
		return(writer);
	}

    /* Start the document with the xml default for the version,
	* encoding ISO 8859-1 and the default for the standalone
    * declaration. */
	rc = xmlTextWriterStartDocument(writer, NULL, MY_ENCODING, NULL);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterStartDocument\n");
		return(writer);
	}

	/* Start an element. */
	rc = xmlTextWriterStartElementNS(writer,
					 BAD_CAST NULL, BAD_CAST "quakeml", BAD_CAST "http://quakeml.ethz.ch/ns/eventlist");
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterStartElementNS\n");
		return(writer);
	}


	return(writer);

}



/**
 * nll2qml:
 * @uri: the output URI
 *
 * test the xmlWriter interface when writing to a new file
 */
void startQuakeMLEventElement(xmlTextWriterPtr writer, const char* uniqueId)
{
	int rc;


	/* Start an element. */
	rc = xmlTextWriterStartElementNS(writer,
					 NULL, BAD_CAST "event", BAD_CAST "http://quakeml.ethz.ch/ns/eventlist/event");
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterStartElementNS\n");
		return;
	}

	/* Add an attribute. */
	rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "unique_id", BAD_CAST uniqueId);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteAttribute\n");
		return;
	}


}



/**
 * writeQuakeMLLocationNLL:
 *
 */
void writeQuakeMLLocationNLL(xmlTextWriterPtr writer,
			  const char* uniqueId, const char* timezone,
			  HypoDesc* phypo, ArrivalDesc* parrivals, int narrivals,
			  int iWriteArrivals, int iWriteMinimal,
			  GridDesc* pgrid, int n_proj)
{
	int rc;
	xmlChar *tmp = NULL;


	/* Start an element. */
	rc = xmlTextWriterStartElementNS(writer,
					 NULL, BAD_CAST "location", BAD_CAST "http://quakeml.ethz.ch/ns/eventlist/location");
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterStartElementNS\n");
		return;
	}
	/* Add an attribute. */
	rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "main", BAD_CAST "true");
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteAttribute\n");
		return;
	}
	/* Add an attribute. */
	rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "unique_id", BAD_CAST uniqueId);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteAttribute\n");
		return;
	}

	char comment[FILENAME_MAX];
	strcat(comment, "NLL file root: ");
	strcat(comment, phypo->fileroot);
	strcat(comment, " ");
	/* Write a comment as child of current element.
	* Please observe, that the input to the xmlTextWriter functions
	* HAS to be in UTF-8, even if the output XML is encoded in iso-8859-1 */
	tmp = ConvertInput(comment, MY_ENCODING);
	rc = xmlTextWriterWriteComment(writer, tmp);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteComment\n");
		return;
	}
	if (tmp != NULL) xmlFree(tmp);


    	// write QuakeML location chidren
	writeQuakeMLOriginTime(writer, timezone, phypo->year, phypo->month, phypo->day, phypo->hour, phypo->min, phypo->sec);
	writeQuakeMLLatitude(writer, phypo->dlat, "degrees", 0.0);
	writeQuakeMLLongitude(writer, phypo->dlong, "degrees", 0.0);
	writeQuakeMLDepth(writer, phypo->depth, "km", 0.0);
	writeQuakeMLMagnitude(writer, phypo->amp_mag, "M", 0.0);
	writeQuakeMLMagnitude(writer, phypo->dur_mag, "MD", 0.0);
	writeQuakeMLRegion(writer, "The World");




}



/**
 * writeQuakeMLOriginTime:
 *
 */
void writeQuakeMLOriginTime(xmlTextWriterPtr writer, const char* timezone,
			    int year, int month, int day, int hour, int min, double sec)
{
	int rc;


	/* Start an element. */
	rc = xmlTextWriterStartElement(writer, BAD_CAST "origin-time");
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterStartElement\n");
		return;
	}
	/* Add an attribute. */
	rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "timezone", BAD_CAST timezone);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteAttribute\n");
		return;
	}

	/* Write an element as child. */
	rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "year", "%d", year);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteFormatElement\n");
		return;
	}

	/* Write an element as child. */
	rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "month", "%d", month);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteFormatElement\n");
		return;
	}

	/* Write an element as child. */
	rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "day", "%d", day);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteFormatElement\n");
		return;
	}

	/* Write an element as child. */
	rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "hour", "%d", hour);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteFormatElement\n");
		return;
	}

	/* Write an element as child. */
	rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "min", "%d", min);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteFormatElement\n");
		return;
	}

	/* Write an element as child. */
	rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "sec", "%lf", sec);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteFormatElement\n");
		return;
	}


	/* Close the element. */
	rc = xmlTextWriterEndElement(writer);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterEndElement\n");
		return;
	}

}



/**
 * writeQuakeMLLatitude:
 *
 */
void writeQuakeMLLatitude(xmlTextWriterPtr writer, double latitude, const char* units, double error)
{
	int rc;


	/* Start an element. */
	rc = xmlTextWriterStartElement(writer, BAD_CAST "latitude");
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterStartElement\n");
		return;
	}
	/* Add an attribute. */
	rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "unit", BAD_CAST units);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteAttribute\n");
		return;
	}
	/* Add an attribute. */
	rc = xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "error", "%lf", error);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteFormatAttribute\n");
		return;
	}

	/* Write formatted CDATA. */
	rc = xmlTextWriterWriteFormatRaw(writer, "%lf", latitude);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteFormatRaw\n");
		return;
	}

	/* Close the element. */
	rc = xmlTextWriterEndElement(writer);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterEndElement\n");
		return;
	}

}



/**
 * writeQuakeMLLongitude:
 *
 */
void writeQuakeMLLongitude(xmlTextWriterPtr writer, double longitude, const char* units, double error)
{
	int rc;


	/* Start an element. */
	rc = xmlTextWriterStartElement(writer, BAD_CAST "longitude");
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterStartElement\n");
		return;
	}
	/* Add an attribute. */
	rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "unit", BAD_CAST units);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteAttribute\n");
		return;
	}
	/* Add an attribute. */
	rc = xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "error", "%lf", error);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteFormatAttribute\n");
		return;
	}

	/* Write formatted CDATA. */
	rc = xmlTextWriterWriteFormatRaw(writer, "%lf", longitude);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteFormatRaw\n");
		return;
	}

	/* Close the element. */
	rc = xmlTextWriterEndElement(writer);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterEndElement\n");
		return;
	}

}



/**
 * writeQuakeMLDepth:
 *
 */
void writeQuakeMLDepth(xmlTextWriterPtr writer, double depth, const char* units, double error)
{
	int rc;


	/* Start an element. */
	rc = xmlTextWriterStartElement(writer, BAD_CAST "depth");
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterStartElement\n");
		return;
	}
	/* Add an attribute. */
	rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "unit", BAD_CAST units);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteAttribute\n");
		return;
	}
	/* Add an attribute. */
	rc = xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "error", "%lf", error);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteFormatAttribute\n");
		return;
	}

	/* Write formatted CDATA. */
	rc = xmlTextWriterWriteFormatRaw(writer, "%lf", depth);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteFormatRaw\n");
		return;
	}

	/* Close the element. */
	rc = xmlTextWriterEndElement(writer);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterEndElement\n");
		return;
	}

}



/**
 * writeQuakeMLMagnitude:
 *
 */
void writeQuakeMLMagnitude(xmlTextWriterPtr writer, double magnitude, const char* units, double error)
{
	int rc;


	/* Start an element. */
	rc = xmlTextWriterStartElement(writer, BAD_CAST "magnitude");
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterStartElement\n");
		return;
	}
	/* Add an attribute. */
	rc = xmlTextWriterWriteAttribute(writer, BAD_CAST "unit", BAD_CAST units);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteAttribute\n");
		return;
	}
	/* Add an attribute. */
	rc = xmlTextWriterWriteFormatAttribute(writer, BAD_CAST "error", "%lf", error);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteFormatAttribute\n");
		return;
	}

	/* Write formatted CDATA. */
	rc = xmlTextWriterWriteFormatRaw(writer, "%lf", magnitude);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteFormatRaw\n");
		return;
	}

	/* Close the element. */
	rc = xmlTextWriterEndElement(writer);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterEndElement\n");
		return;
	}

}



/**
 * writeQuakeMLRegion:
 *
 */
void writeQuakeMLRegion(xmlTextWriterPtr writer, const char* region)
{
	int rc;


	/* Write an element as child. */
	rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "region", "%s", region);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterWriteFormatElement\n");
		return;
	}


}



/**
 * endQuakeMLDocument:
 *
 */
void endQuakeMLDocument(xmlTextWriterPtr writer)
{
	int rc;

    /* Here we could close remaining open elements using the
	* function xmlTextWriterEndElement, but since we do not want to
	* write any other elements, we simply call xmlTextWriterEndDocument,
    * which will do all the work. */
	rc = xmlTextWriterEndDocument(writer);
	if (rc < 0) {
		printf("nll2qml: Error at xmlTextWriterEndDocument\n");
		return;
	}

	xmlFreeTextWriter(writer);

}






/* ----------------------------------------------------------------------------------- */



/**
 * this initialize the library and check potential ABI mismatches
 * between the version it was compiled for and the actual shared
 * library used.
 */
int
		initLibXml(void)
{

	LIBXML_TEST_VERSION

			return 0;

}


/**
 * Cleanup function for the XML library.
 */
int
		cleanUpLibXml(void)
{

    /*
	* Cleanup function for the XML library.
    */
	xmlCleanupParser();

	return 0;
}


/**
 * ConvertInput:
 * @in: string in a given encoding
 * @encoding: the encoding used
 *
 * Converts @in into UTF-8 for processing with libxml2 APIs
 *
 * Returns the converted UTF-8 string, or NULL in case of error.
 */
xmlChar *
		ConvertInput(const char *in, const char *encoding)
{
	xmlChar *out;
	int ret;
	int size;
	int out_size;
	int temp;
	xmlCharEncodingHandlerPtr handler;

	if (in == 0)
		return 0;

	handler = xmlFindCharEncodingHandler(encoding);

	if (!handler) {
		printf("ConvertInput: no encoding handler found for '%s'\n",
		       encoding ? encoding : "");
		return 0;
	}

	size = (int) strlen(in) + 1;
	out_size = size * 2 - 1;
	out = (unsigned char *) xmlMalloc((size_t) out_size);

	if (out != 0) {
		temp = size - 1;
		ret = handler->input(out, &out_size, (const xmlChar *) in, &temp);
		if ((ret < 0) || (temp - size + 1)) {
			if (ret < 0) {
				printf("ConvertInput: conversion wasn't successful.\n");
			} else {
				printf
						("ConvertInput: conversion wasn't successful. converted: %i octets.\n",
						temp);
			}

			xmlFree(out);
			out = 0;
		} else {
			out = (unsigned char *) xmlRealloc(out, out_size + 1);
			out[out_size] = 0;  /*null terminating out */
		}
	} else {
		printf("ConvertInput: no mem\n");
	}

	return out;
}



#else
int main(void) {
	fprintf(stderr, "ERROR - libxml writer or output support not compiled in.\n");
	exit(1);
}
#endif

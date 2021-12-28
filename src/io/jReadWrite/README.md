# jReadWrite
JSON Library for use with [Mbed-OS](https://os.mbed.com/)

Credit to [Tony Wilk](https://www.codeproject.com/Members/tonywilk)

## Introduction
jReadWrite is an Mbed library composed of [jWrite]( https://github.com/jonaskgandersson/jWrite) and [jRead](https://github.com/jonaskgandersson/jRead).

**jWrite** is a simple way of writing JSON to a char buffer in C, directly from native variables. It manages the output buffer so you don't overrun, it handles all the fiddly quotes, brackets and commas and reports where you have tried to create invalid JSON.

**jRead** is a simple in-place JSON element reader, it maintains the input JSON as unaltered text and allows queries to be made on it directly. Written in C, jRead is small, fast and does not allocate any memory - making it an ideal choice for embedded applications or those cases where all you want is to read a few values from a chunk of JSON without having to learn a large API and have your program generate lots of interlinked structures to represent the JSON.

## Getting Started
Here is an Mbed example project [mbed-json-example](https://github.com/jonaskgandersson/mbed-json-example)
#### Adding a library to your project
Use `mbed add` to add the latest revision of a library:

```
$ mbed add https://github.com/jonaskgandersson/jReadWrite
```

Use the `URL#hash` format to add a library from a URL at a specific revision hash:

```
$ mbed add https://github.com/jonaskgandersson/jReadWrite/#6dd542b2814768fdc6ca0aaea7529dfcb9bee37f
```
### jWrite
Jumping straight in:
```c
struct jWriteControl jwc;
jwOpen( &jwc, buffer, buflen, JW_OBJECT, JW_PRETTY );  // open root node as object
jwObj_string( &jwc, "key", "value" );                  // writes "key":"value"
jwObj_int( &jwc, "int", 1 );                           // writes "int":1
jwObj_array( &jwc, "anArray");                         // start "anArray": [...] 
    jwArr_int( &jwc, 0 );                              // add a few integers to the array
    jwArr_int( &jwc, 1 );
    jwArr_int( &jwc, 2 );
jwEnd( &jwc );                                         // end the array
err= jwClose( &jwc );                                  // close root object - done
```

which results in:
```c
{
    "key": "value",
    "int": 1,
    "anArray": [
        0,
        1,
        2
    ]
}
```
The output is prettyfied (it's an option) and has all the { } [ ] , : " characters in the right place.

#### To Be, or Not To Be, GLOBAL
For many applications it is a lot simpler to have one global (static) instance of a structure which can be used for jWrite, it makes the API calls easy to type in - you don't have to supply a reference every time.

The above example with JW_GLOBAL_CONTROL_STRUCT set to true in mbed_lib.json:
```c
jwOpen( buffer, buflen, JW_OBJECT, JW_PRETTY );  // open root node as object
jwObj_string( "key", "value" );                  // writes "key":"value"
jwObj_int( "int", 1 );                           // writes "int":1
jwObj_array( "anArray");                         // start "anArray": [...] 
    jwArr_int( 0 );                              // add a few integers to the array
    jwArr_int( 1 );
    jwArr_int( 2 );
jwEnd();                                         // end the array
err= jwClose();                                  // close root object - done
```
However, that is not very flexible and does not allow for multiple uses of jWrite functions at the same time.

### jRead
Basically, the API is a single function call:

```c
jRead( pJsonText, pQueryString, pOutputStruct );
```
The Query String defines how to traverse the JSON to end up with a value which is returned to pOutputStruct . Elements of a query can be as follows:

```c
"{'keyname'"           Object element "keyname", returns value of that key
"{NUMBER"              Object element[NUMBER], returns keyname of that element
"[NUMBER"              Array element[NUMBER], returns value from array
```

Examples of What jRead Does
Assume we have a simple JSON string we have read into some buffer pJson like:

```c
{
    "astring":"This is a string",
    "anumber":42,
    "myarray":[ "zero", 1, {"description":"element 2"}, null ],
    "yesno":true,
    "PI":"3.1415926",
    "foo":null
}
```
we can get the value of any element using jRead, e.g.,

```c
struct jReadElement result;
jRead( pJson, "{'astring'", &result );
```
In this case, the result would be:

```c
result.dataType= JREAD_STRING
result.elements= 1
result.byteLen=  16
result.pValue -> "This is a string"
```
Note that nothing is copied, the result structure simply gives you a pointer to the start of the element and its length.

The helper functions make it easy to extract single values into C variables like:

```c
jRead_string( pJson, "{'astring'", destString, MAXLEN, NULL );	// This is a string
my_int = jRead_int( pJson, "{'myarray'[1", NULL );		// 1
my_long = jRead_long( pJson, "{'anumber'", NULL );		// 42
my_double = jRead_double( pJson, "{'PI'", NULL );		// 3.1415926
```
Note that the 'int' helper does some 'type coersion' - since everything is a string anyway, calling jReadInt() will always return a value: it will return 42 for the JSON 42 or "42" and returns zero for "foo". It also returns 1 for true and 0 for false or null.

## Compiled with

* Mbed CLI [GNU Arm Embedded toolchain (GCC)](https://developer.arm.com/open-source/gnu-toolchain/gnu-rm/downloads)

## Contributing

When contributing to this repository, please first discuss the change you wish to make via issue, email, or any other method with the owners of this repository before making a change.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/jonaskgandersson/jReadWrite/tags). 

## Authors

* **Jonas Andersson** - *Initial work* - [Tony Wilk](https://www.codeproject.com/Members/tonywilk)


## License

This project is licensed under the Apache License - see the [LICENSE](https://github.com/jonaskgandersson/jReadWrite/blob/master/LICENSE) file for details

## Acknowledgments

* Hat tip to [Tony Wilk](https://www.codeproject.com/Members/tonywilk) for sharing his great code.

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2009-2010 Michael Specht.
# module is distributed under the terms of the GNU General Public License

'''
parse string into token list
parse token list
'''

"""
This file is part of p3d.

    p3d is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    p3d is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
"""


import bisect, copy, re
from collections import defaultdict as ddict

class Parser:
    # atomic tokens are retained during tokenization
    ATOMIC_TOKENS = sorted([
        '!=', '<>', '..', '||', '&&', '<=', '>=', '==', '=', 
        '<', '>', '!', '&', '|', '^', ',', '(', ')', '{', '}', '-', ':'
        ], key=lambda x : len(x), reverse = True)
    
    LOGICAL_OPERATORS = {
        'not': ['not', '!'],
        'and': ['and', '&&', '&'],
        'xor': ['xor', '^'],
        'or': ['or', '||', '|'],
        'minus': ['-']
    }
    
    LOGICAL_OPERATOR_PRECEDENCE = ['and', 'xor', 'or', 'minus']
    
    ARITHMETIC_OPERATORS = {
        'e': [['is'], ['=='], ['=']],
        'ne': [['is', 'not'], ['!='], ['<>']],
        'l': [['<']],
        'le': [['<=']],
        'g': [['>']],
        'ge': [['>=']]
    }
    
    VALUE_SEPARATOR = [',']
    
    VALUE_RANGE = ['to', '..']
    
    def __init__(self, repository, meta = None):
        if meta == None: meta = {}
        if not 'aliases' in meta: meta['aliases'] = {}
        if not 'functions' in meta: meta['functions'] = {}
        if not 'caseSensitive' in meta: meta['caseSensitive'] = set([])
        self.repository = repository
        self.initRepository(meta['aliases'], meta['functions'], meta['caseSensitive'])
        
    def initRepository(self, aliases, functions, caseSensitive):
        # convert all repository hash keys
        # to lowercase because if something comes in as ALA, it will never ever
        # be recognized by this parser which lowercases all input first
        repositoryCopy = dict()
        for key in self.repository:
            lowerKey = key.lower()
            repositoryCopy[lowerKey] = dict()
            for subKey in self.repository[key]:
                changedSubKey = subKey
                if (type(subKey) == str) and (not (key in caseSensitive)):
                    changedSubKey = changedSubKey.lower()
                repositoryCopy[lowerKey][changedSubKey] = self.repository[key][subKey]
        self.repository = repositoryCopy
        
        self.repositoryMeta = dict()
        # merge all sets into the global 'all' set
        self.repositoryMeta['globalAll'] = set()
        for values in self.repository.values():
            for subset in values.values():
                self.repositoryMeta['globalAll'] |= subset
                
        # create the 'local all' sets if they differ from the global all set
        self.repositoryMeta['localAll'] = dict()
        for key in self.repository.keys():
            localAll = set()
            for subset in self.repository[key].values():
                localAll |= subset
            if (localAll == self.repositoryMeta['globalAll']):
                localAll = self.repositoryMeta['globalAll']
            self.repositoryMeta['localAll'][key] = localAll
        
        # index key/value sets
        self.repositoryMeta['index'] = dict()
        for key, value in self.repository.items():
            self.repositoryMeta['index'][key] = sorted(list(value.keys()))
        
        # determine key types for kv sets
        self.repositoryMeta['type'] = dict()
        for key, value in self.repository.items():
            keyType = None
            for subKey in value.keys():
                if (keyType == None):
                    keyType = type(subKey)
                else:
                    if (keyType != type(subKey)):
                        raise Exception("Inconsistent types accross the " + key + "hash.")
            self.repositoryMeta['type'][key] = keyType
            ''' e.g. non-aa-resname = <class 'str'> '''
            
        # add repository keys to alias hash
        for key in self.repository.keys():
            aliases[key] = key
        aliases['all'] = 'all'
        self.repositoryMeta['type']['all'] = str
        self.repositoryMeta['aliases'] = aliases
        
        # promote aliases to array of ([token list], key) tuples (so that 'atom type' is possible)
        aliasList = []
        for key in self.repositoryMeta['aliases']:
            keySplit = self.tokenize(key.lower())
            aliasList.append((keySplit, self.repositoryMeta['aliases'][key]))
        self.repositoryMeta['aliases'] = aliasList  
        
        # parse function definitions: 'length in {items: set} is {length: int}'
        self.repositoryMeta['functions'] = []
        for pattern in functions:
            tokens = self.tokenize(pattern)
            description = []
            while len(tokens) > 0:
                if tokens[0] == '{':
                    tokens.pop(0)
                    argumentId = tokens.pop(0)
                    tokens.pop(0)
                    argumentType = tokens.pop(0)
                    if argumentType == 'float':
                        argumentType = float
                    elif argumentType == 'int':
                        argumentType = int
                    elif argumentType == 'set':
                        argumentType = set
                    elif argumentType == 'value':
                        # value of key
                        tokens.pop(0)
                        valueKey = tokens.pop(0)
                        if not (valueKey in self.repository):
                            raise Exception("Invalid argument type value of " + valueKey + ".")
                        argumentType = valueKey
                    else:
                        raise Exception("Invalid argument type " + argumentType + ".")
                    tokens.pop(0)
                    description.append((argumentId, argumentType))
                else:
                    description.append(tokens.pop(0).lower())
            self.repositoryMeta['functions'].append({'pattern': description, 'function': functions[pattern]})
            self.repositoryMeta['caseSensitive'] = caseSensitive
        return
        
    def tokenize(self, term):
        replacementTable = {}
        # replace each atomic token with `1 plus whitespace etc.
        for index, token in enumerate(self.ATOMIC_TOKENS):
            # don't replace minus (-) if it is followed by a digit (0-9) or a decimal point (.)
            tokenRegex = None
            if token == '-':
                tokenRegex = '\-[^0-9\.]'
            placeHolder = "`{0}".format(index)
            if tokenRegex:
                # use a regex for minus sign
                term = re.sub(tokenRegex, ' ' + placeHolder + ' ', term)
            else:
                # use simple search & replace
                term = term.replace(token, ' ' + placeHolder + ' ')
            replacementTable[placeHolder] = token
        termSplit = term.split()
        
        # replace `1 etc. with their original meaning
        for index, token in enumerate(termSplit):
            if token in replacementTable.keys():
                termSplit[index] = replacementTable[token]
                
        return termSplit
        
    def lowerList(self, inlist):
        '''
        lowercase all entries in list
        '''
        llist = []
        for item in inlist:
            llist.append(item.lower())
        return llist
    
    def parseSetKey(self, tokens):
        '''
        see if we find any set key. if something is found, remember
        the token count of this operator in maxAliasLength (i.e. ATOM TYPE has two
        tokens, while ATYPE has only one), then pick the one with the highest token
        count. this is to ensure that the longest possible key is picked.
        '''
        localTokens = copy.copy(tokens)
        resultKey = None
        maxAliasLength = 0
        
        for tokensAndPair in self.repositoryMeta['aliases']:
            keyAlias = tokensAndPair[0]
            tokenSlice = self.lowerList(localTokens[0:len(keyAlias)])
            if (tokenSlice == keyAlias) and (resultKey == None or (len(keyAlias) > maxAliasLength)):
                # we have found something!
                resultKey = tokensAndPair[1]
                maxAliasLength = len(keyAlias)
                    
        if resultKey == None:
            localTokens = tokens
        else:
            del localTokens[0:maxAliasLength]
            
        return resultKey, localTokens
    
    def parseLogicalOperator(self, tokens):
        '''
        see if we find any logical operator. 
        '''
        if len(tokens) == 0:
            return None, tokens
            
        localTokens = copy.copy(tokens)
        
        for key in self.LOGICAL_OPERATORS:
            for keyAlias in self.LOGICAL_OPERATORS[key]:
                tokenSlice = localTokens[0].lower()
                if tokenSlice == keyAlias:
                    # we have found something!
                    localTokens.pop(0)
                    return key, localTokens
                    
        return None, tokens
    
    def parseArithmeticOperator(self, tokens):
        '''
        see if we find any arithmetic operator. if something is found, remember
        the token count of this operator in maxAliasLength (i.e. IS NOT has two
        tokens, while IS has only one), then pick the one with the highest token
        count. this is to ensure that IS NOT will not be interpreted as IS, followed
        by a NOT.
        '''
        localTokens = copy.copy(tokens)
        resultKey = None
        maxAliasLength = 0
        
        for key in self.ARITHMETIC_OPERATORS:
            for keyAlias in self.ARITHMETIC_OPERATORS[key]:
                tokenSlice = self.lowerList(localTokens[0:len(keyAlias)])
                if (tokenSlice == keyAlias) and (resultKey == None or (len(keyAlias) > maxAliasLength)):
                    # we have found something!
                    resultKey = key
                    maxAliasLength = len(keyAlias)
                    
        if resultKey == None:
            localTokens = tokens
        else:
            del localTokens[0:maxAliasLength]
        
        return resultKey, localTokens
    
    def parseSingleValue(self, tokens, valueType, key):
        localTokens = copy.copy(tokens)
        
        # make sure that this is not a set key
        result, localTokens = self.parseSetKey(localTokens)
        if result != None:
            return None, tokens
        
        # make sure that this is not a logical operator
        result, localTokens = self.parseLogicalOperator(localTokens)
        if result != None:
            return None, tokens
        
        leadingMinus = False
        if (valueType == int or valueType == float):
            if (localTokens[0] == '-'):
                leadingMinus = True
                localTokens.pop(0)
        value = None
        try:
            value = valueType(localTokens.pop(0))
        except:
            return None, tokens
            
        if (valueType == str) and (not (key in self.repositoryMeta['caseSensitive'])): 
            value = value.lower()
            
        if leadingMinus:
            value = -value
        return value, localTokens
    
    def parseValue(self, tokens, valueType, key):
        '''
        see if we can parse a value. A value may be described in the form of
        1 .. 4 , 6 , 10 - 11 , 14 to 18
        the return value is an array of arrays:
        [[1, 4], 6, [10, 11], [14, 18]]
        with a single value 4, the array would be: [[4]]
        '''
        localTokens = copy.copy(tokens)
        
        valueList = []
        
        while len(localTokens) > 0:
            valueStart, localTokens = self.parseSingleValue(localTokens, valueType, key)
            if valueStart == None:
                break
                
            valueEnd = None
            if len(localTokens) > 1 and (localTokens[0].lower() in self.VALUE_RANGE):
                localTokens.pop(0)
                valueEnd, localTokens = self.parseSingleValue(localTokens, valueType, key)
                if valueEnd == None:
                    raise Exception("Range operator but no upper value defined.")
                
            thisValue = [valueStart]
            if valueEnd != None:
                thisValue.append(valueEnd)
            valueList.append(thisValue)
            
            if not (len(localTokens) > 1 and localTokens[0].lower() in self.VALUE_SEPARATOR):
                break
            localTokens.pop(0)
        
        return valueList, localTokens
    
    def parseSingleFunction(self, tokens, function):
        localTokens = copy.copy(tokens)
        localPattern = copy.copy(function['pattern'])
        arguments = {}
        while len(localPattern) > 0 and len(localTokens) > 0:
            part = localPattern.pop(0)
            if (type(part) == str):
                if part != localTokens.pop(0).lower():
                    break
            else:
                argumentId = part[0]
                argumentType = part[1]
                argumentValue = None
                if (argumentType == set):
                    argumentValue, localTokens = self.parseOneSet(localTokens, parseFunction = False)
                elif (type(argumentType) == str):
                    # it's a value of a certain type!
                    argumentValue, localTokens = self.parseSingleValue(localTokens, str, argumentType)
                else:
                    argumentValue, localTokens = self.parseSingleValue(localTokens, argumentType, None)
                if argumentValue == None:
                    break
                arguments[argumentId] = argumentValue
            if len(localPattern) == 0:
                result = function['function'](**arguments)
                return result, localTokens
                
        return None, tokens
    
    def parseFunction(self, tokens):
        localTokens = copy.copy(tokens)
        for function in self.repositoryMeta['functions']:
            result, localTokens = self.parseSingleFunction(localTokens, function)
            if result != None:
                return result, localTokens
        return None, tokens
        
    def parseOneSet(self, tokens, parseFunction = True):
        '''
        try to parse the following:
          'not', SET
        | '(', SET, ')'
        | {0}
        | SPECIAL CONSTRUCT (like WITHIN {length} OF {items})
        | SET_KEY [[, AOP], VALUE]
        
        this function really expects a set, not empty string is allowed here
        return result set or None, tokens
        '''
        localTokens = copy.copy(tokens)
        
        if (len(localTokens) == 0):
            return None, localTokens
            
        # try: not set
        if (localTokens[0].lower() in self.LOGICAL_OPERATORS['not']):
            localTokens.pop(0)
            result, localTokens = self.parseOneSet(localTokens)
            result = self.repositoryMeta['globalAll'] - result
            return result, localTokens
            
        # try (set)
        if (localTokens[0] == '('):
            localTokens.pop(0)
            depth = 1
            subTokens = []
            while (depth > 0):
                subToken = localTokens.pop(0)
                if (subToken == '('):
                    depth += 1
                elif (subToken == ')'):
                    depth -= 1
                if not ((subToken == ')') and (depth == 0)):
                    subTokens.append(subToken)
            result, subTokens = self.parseSet(subTokens)
            return result, localTokens
            
        # try {0}
        if (localTokens[0] == '{'):
            localTokens.pop(0)
            id = int(localTokens.pop(0))
            localTokens.pop(0)
            return self.parameters[id], localTokens
            
        if parseFunction:
            # try functions
            result, localTokens = self.parseFunction(localTokens)
            if (result != None):
                return result, localTokens
        
        # try set key (with operator and values?)
        key, localTokens = self.parseSetKey(localTokens)
        if key != None:
            keyType = self.repositoryMeta['type'][key]
            operator, localTokens = self.parseArithmeticOperator(localTokens)
            value, localTokens = self.parseValue(localTokens, keyType, key)
            if (value == None and operator != None):
                raise Exception("Operator but no value defined.")
            result = self.evaluate(key, operator, value)
            return result, localTokens
            
        return None, tokens
            
    def parseSet(self, tokens):
        if len(tokens) == 0:
            return tokens, set()
            
        localTokens = copy.copy(tokens)
            
        # terms includes sets and operators
        terms = []
        
        while len(localTokens) > 0:
            localSet, localTokens = self.parseOneSet(localTokens)
            if localSet == None:
                break
                
            terms.append(localSet)
            
            if len(localTokens) == 0:
                break
                
            logicalOperator, localTokens = self.parseLogicalOperator(localTokens)
            if logicalOperator == None:
                # raise Exception("Invalid logical operator: " + localTokens[0] + ".")
                # instead of raising an exception, we just silently assume 'and' here
                # example: 'oxygen within 2 of backbone' probably means 'oxygen 
                # AND within 2 of backbone' (we want all atoms that are oxygens 
                # and within 2 Angstrom of the backbone
                logicalOperator = 'and'
            terms.append(logicalOperator)
            
        # now figure out operator precedence and combine collected sets
        for operator in self.LOGICAL_OPERATOR_PRECEDENCE:
            # search for operator in term list
            while (1):
                i = 1
                while (i < len(terms) and terms[i] != operator):
                    i += 2
                if (i < len(terms)):
                    # found an operator!
                    result = None
                    leftSet = terms[i - 1]
                    rightSet = terms[i + 1]
                    
                    if operator == 'and':
                        result = leftSet & rightSet
                    elif operator == 'xor':
                        result = leftSet ^ rightSet
                    elif operator == 'or':
                        result = leftSet | rightSet
                    elif operator == 'minus':
                        result = leftSet - rightSet

                    # delete three items and replace by result
                    del terms[(i-1):(i+2)]
                    terms.insert(i-1, result)
                else:
                    break
                    
        return terms[0], localTokens
    
    def evaluate(self, key, operator = None, value = None, recursion = False):
        if key == 'all':
            return self.repositoryMeta['globalAll']
        
        if not (key in self.repository):
            raise Exception("Unknown identifier '" + key + "'.")
            
        if (value == None or len(value) == 0):
            if (operator == None):
                return self.repositoryMeta['localAll'][key]
            else:
                raise Exception("Operator but no value defined.")
            
        if (operator == None):
            operator = 'e'

        isComplexValue = (len(value) > 1) or (len(value) > 0 and len(value[0]) > 1)
        if (not (operator == 'e' or operator == 'ne')) and isComplexValue:
            raise Exception("Value ranges are only allowed with = or !=")
            
        if (operator == 'e' or operator == 'ne'):
            # handle = and !=, value may be a range
            result = set()
            
            for valueSingleOrRange in value:
                if len(valueSingleOrRange) == 1:
                    # single value, pull out items
                    if (valueSingleOrRange[0] in self.repository[key]):
                        result |= self.repository[key][valueSingleOrRange[0]]
                else:
                    # we have a range, recursively call this evaluate
                    subResult = self.evaluate(key, 'ge', [[valueSingleOrRange[0]]])
                    subResult &= self.evaluate(key, 'le', [[valueSingleOrRange[1]]])
                    result |= subResult
                        
            # negate if necessary
            if (operator == 'ne'):
                result = self.repositoryMeta['localAll'][key] - result
                    
            return result
        else:
            # don't handle range values, this is l, le, g, or ge
            value = value[0][0]
            index = bisect.bisect_left(self.repositoryMeta['index'][key], value)
            indexAtValue = False
            if (index < len(self.repositoryMeta['index'][key])):
                indexAtValue = (self.repositoryMeta['index'][key][index] == value)
            if (not indexAtValue) and (operator == 'l' or operator == 'le'):
                index -= 1
            # shift index if > or <
            if (indexAtValue):
                if (operator == 'l'):
                    index -= 1
                    indexAtValue = False
                elif (operator == 'g'):
                    index += 1
                    indexAtValue = False

            result = set()
            
            valueRange = None
            if (operator == 'g' or operator == 'ge'):
                valueRange = range(index, len(self.repositoryMeta['localAll'][key]))
            else:
                valueRange = range(0, index + 1)
                
            for step in valueRange:
                if (step >= 0 and step < len(self.repositoryMeta['index'][key])):
                    result |= self.repository[key][self.repositoryMeta['index'][key][step]]
                    
            return result
    
    def parse(self, *args):
        term = ''
        parameters = []
        for item in args:
            if type(item) == str: 
                term += ' ' + item + ' '
            else:
                if (type(item) != set): 
                    if type(item) == list:
                        item = set(item)
                    else:
                        item = set([item])
                    
                term += ' {' + str(len(parameters)) + '} '
                parameters.append(item)
        
        self.parameters = parameters
        tokens = self.tokenize(term)
        result, tokens = self.parseSet(tokens)
        # tokens should be empty now
        return result
    
def strlen_e(items=set(), length=0.0):
    result = set()
    for item in items:
        if len(item) == length:
            result.add(item)
    return result

def test():
    repository = dict()
    repository['letter'] = dict()
    repository['letter'][''] = set(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 
                                            'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 
                                            's', 't', 'u', 'v', 'w', 'x', 'y', 'z'])
    repository['digit'] = dict()
    repository['digit'][''] = set(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])
    repository['vowel'] = dict()
    repository['vowel'][''] = set(['a', 'e', 'i', 'o', 'u'])
    repository['even'] = dict()
    repository['even'][''] = set(['0', '2', '4', '6', '8'])
    repository['odd'] = dict()
    repository['odd'][''] = set(['1', '3', '5', '7', '9'])
    repository['mod3'] = dict()
    repository['mod3'][0] = set(['0', '3', '6', '9'])
    repository['mod3'][1] = set(['1', '4', '7'])
    repository['mod3'][2] = set(['2', '5', '8'])
    repository['prime'] = dict()
    repository['prime'][2] = set(['2'])
    repository['prime'][3] = set(['3'])
    repository['prime'][5] = set(['5'])
    repository['prime'][7] = set(['7'])
    repository['prime'][11] = set(['11'])
    repository['prime'][13] = set(['13'])
    repository['prime'][17] = set(['17'])
    
    aliases = dict()
    aliases['prime number'] = 'prime'
    
    functions = dict()
    functions['length is {length: float} in {items: set}'] = strlen_e
    
    meta = dict()
    meta['aliases'] = aliases
    meta['functions'] = functions
    parser = Parser(repository, meta)
    
    queries = [
        'letter and vowel',
        'letter vowel',
        #'prime is not 5',
        #'Prime >= 4',
        #'PRIME >= 4',
        #'prime number >= 4',
        #'(length is 1 in prime) or vowel',
        #'length is 1 in vowel or prime',
        #'length is 1 in (prime or vowel)',
        #'prime is 3..4,5,8..16 - digit',
    ]
    for query in queries:
        #try:
            print(query + ':', parser.parse(query))
        #except Exception as e:
            #print(query + ': Error:', e)
    '''
    print(parser.parse("(())"))
    print(parser.parse("(all)"))
    print(parser.parse("digit and not even"))
    print(parser.parse("(digit and not even) or vowel"))
    print(parser.parse("digit and (not even or vowel)"))
    print(parser.parse("(digit and not even) or not vowel"))
    print(parser.parse("digit and not mod3 is 0"))
    print(parser.parse("not mod3 is 0"))
    print(parser.parse("mod3 is not 0"))
    print(parser.parse("digit and even or vowel"))
    print(parser.parse("vowel or even and digit"))
    print(parser.parse("(vowel or even) and digit"))
    print(parser.parse('prime == 5'))
    print(parser.parse('prime != 5'))
    print(parser.parse('prime > 5'))
    print(parser.parse('prime >= 5'))
    print(parser.parse('prime < 5'))
    print(parser.parse('prime <= 5'))
    print(parser.parse('prime == 4'))
    print(parser.parse('prime != 4'))
    print(parser.parse('prime > 4'))
    print(parser.parse('prime >= 4'))
    print(parser.parse('prime < 4'))
    print(parser.parse('prime <= 4'))
    print(parser.parse('prime > 100'))
    print(parser.parse('prime <= -100'))
    print(parser.parse('prime >= -100'))
    print(parser.parse('prime < 100'))
    print(parser.parse('prime is 3 to 6'))
    print(parser.parse('prime 3..6,11'))
    print(parser.parse('prime is not 3 to 6'))
    print(parser.parse('prime != 3..6,11'))
    print(parser.parse('prime'))
    '''
    exit(0)

        
if __name__ == '__main__':
    print('Running test suite...')
    test()

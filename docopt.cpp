//
//  docopt.cpp
//  docopt
//
//  Created by Jared Grubb on 2013-11-03.
//  Copyright (c) 2013 Jared Grubb. All rights reserved.
//

#include "docopt.h"
#include "docopt_util.h"
#include "docopt_private.h"

#include "docopt_value.h"

#include <vector>
#include <unordered_set>
#include <algorithm>
#include <map>
#include <string>
#include <iostream>
#include <cassert>

using namespace docopt;

DocoptExitHelp::DocoptExitHelp()
: std::runtime_error("Docopt --help argument encountered")
{}

DocoptExitVersion::DocoptExitVersion()
: std::runtime_error("Docopt --version argument encountered")
{}

const char* value::kindAsString(Kind kind)
{
	switch (kind) {
		case Kind::Empty: return "empty";
		case Kind::Bool: return "bool";
		case Kind::Long: return "long";
		case Kind::String: return "string";
		case Kind::StringList: return "string-list";
	}
	return "unknown";
}

void value::throwIfNotKind(Kind expected) const
{
	if (kind == expected)
		return;
	
	std::string error = "Illegal cast to ";
	error += kindAsString(expected);
	error += "; type is actually ";
	error += kindAsString(kind);
	throw std::runtime_error(std::move(error));
}

std::string value::describe() const {
	if (this->isBool()) {
		bool b = this->asBool();
		return b ? "true" : "false";
	} else if (this->isLong()) {
		long v = this->asLong();
        return std::to_string(v);
	} else if (this->isString()) {
		std::string const& str = this->asString();
		return '"' + str + '"';
	} else if (this->isStringList()) {
		auto const& list = this->asStringList();
		std::string result = "[";
		bool first = true;
		for(auto const& el : list) {
            if (! first) {
                result.append(", ");
            }
            first = false;
            result.push_back('"');
            result.append(el);
            result.push_back('"');
		}
		result.push_back(']');
        return result;
	} else {
		return "null";
	}

}

std::ostream& docopt::operator<<(std::ostream& os, value const& val)
{
    os << val.describe();
    return os;
}

#pragma mark -
#pragma mark Pattern types

std::vector<LeafPattern*> Pattern::leaves() {
	std::vector<LeafPattern*> ret;
	collect_leaves(ret);
	return ret;
}

bool Required::match(PatternList& left, std::vector<std::shared_ptr<LeafPattern>>& collected) const
{
	auto l = left;
	auto c = collected;
	
	for(auto const& pattern : fChildren) {
		bool ret = pattern->match(l, c);
		if (!ret) {
			// leave (left, collected) untouched
			return false;
		}
	}
	
	left = std::move(l);
	collected = std::move(c);
	return true;
}

bool LeafPattern::match(PatternList& left, std::vector<std::shared_ptr<LeafPattern>>& collected) const
{
	auto match = single_match(left);
	if (!match.second) {
		return false;
	}
	
	left.erase(left.begin()+match.first);
	
	auto same_name = std::find_if(collected.begin(), collected.end(), [&](std::shared_ptr<LeafPattern> const& p) {
		return p->name()==name();
	});
	if (getValue().isLong()) {
		long val = 1;
		if (same_name == collected.end()) {
			collected.push_back(match.second);
			match.second->setValue(value{val});
		} else if ((**same_name).getValue().isLong()) {
			val += (**same_name).getValue().asLong();
			(**same_name).setValue(value{val});
		} else {
			(**same_name).setValue(value{val});
		}
	} else if (getValue().isStringList()) {
		std::vector<std::string> val;
		if (match.second->getValue().isString()) {
			val.push_back(match.second->getValue().asString());
		} else if (match.second->getValue().isStringList()) {
			val = match.second->getValue().asStringList();
		} else {
			/// cant be!?
		}
		
		if (same_name == collected.end()) {
			collected.push_back(match.second);
			match.second->setValue(value{val});
		} else if ((**same_name).getValue().isStringList()) {
			std::vector<std::string> const& list = (**same_name).getValue().asStringList();
			val.insert(val.begin(), list.begin(), list.end());
			(**same_name).setValue(value{val});
		} else {
			(**same_name).setValue(value{val});
		}
	} else {
		collected.push_back(match.second);
	}
	return true;
}

Option Option::parse(std::string const& options_and_description)
{
    int argCount = 0;
    std::string shortOption, longOption;
    value defaultValue {false};
    
    // Split options_and_description about the double space
    std::string options, description;
    size_t doubleSpace = options_and_description.find("  ");
    if (doubleSpace == std::string::npos)
    {
        // No description
        options = trim(options_and_description);
    }
    else
    {
        options = options_and_description.substr(0, doubleSpace);
        description = trim(options_and_description.substr(doubleSpace));
    }
    
    // Replace equals and commas with spaces
    std::replace(options.begin(), options.end(), ',', ' ');
    std::replace(options.begin(), options.end(), '=', ' ');
    
    // Split options by whitespace
    std::vector<std::string> splitOptions = split(options);
    
    // Walk over them
    for (const std::string &s : splitOptions) {
        if (starts_with(s, "--")) {
            longOption = s;
        } else if (starts_with(s, "-")) {
            shortOption = s;
        } else {
            argCount = 1;
        }
    }
    
    if (argCount) {
        // Extract a default value (e.g. 'foo') from the first substring looks like [default: foo]
        // Unlike Python, we do not require a space after the colon
        const char *prefix = "[default:";
        size_t defaultStartPos = find_case_insensitive(description, prefix);
        if (defaultStartPos != std::string::npos) {
            // Skip past the prefix
            defaultStartPos += strlen(prefix);
            size_t defaultEndPos = description.find("]", defaultStartPos);
            defaultValue = trim(description.substr(defaultStartPos, defaultEndPos - defaultStartPos));
        } else {
            // TODO: generate a warning about missing end ]
        }
    }
	return {std::move(shortOption),
		std::move(longOption),
		argCount,
		std::move(defaultValue)};
}

bool OneOrMore::match(PatternList& left, std::vector<std::shared_ptr<LeafPattern>>& collected) const
{
	assert(fChildren.size() == 1);
	
	auto l = left;
	auto c = collected;
	
	bool matched = true;
	size_t times = 0;
	
	decltype(l) l_;
	bool firstLoop = true;
	
	while (matched) {
		// could it be that something didn't match but changed l or c?
		matched = fChildren[0]->match(l, c);
		
		if (matched)
			++times;
		
		if (firstLoop) {
			firstLoop = false;
		} else if (l == l_) {
			break;
		}
		
		l_ = l;
	}
	
	if (times == 0) {
		return false;
	}
	
	left = std::move(l);
	collected = std::move(c);
	return true;
}

bool Either::match(PatternList& left, std::vector<std::shared_ptr<LeafPattern>>& collected) const
{
	using Outcome = std::pair<PatternList, std::vector<std::shared_ptr<LeafPattern>>>;
	
	std::vector<Outcome> outcomes;
	
	for(auto const& pattern : fChildren) {
		// need a copy so we apply the same one for every iteration
		auto l = left;
		auto c = collected;
		bool matched = pattern->match(l, c);
		if (matched) {
			outcomes.emplace_back(std::move(l), std::move(c));
		}
	}
	
	auto min = std::min_element(outcomes.begin(), outcomes.end(), [](Outcome const& o1, Outcome const& o2) {
		return o1.first.size() < o2.first.size();
	});
	
	if (min == outcomes.end()) {
		// (left, collected) unchanged
		return false;
	}
	
	std::tie(left, collected) = std::move(*min);
	return true;
}

std::pair<size_t, std::shared_ptr<LeafPattern>> Argument::single_match(PatternList const& left) const
{
	std::pair<size_t, std::shared_ptr<LeafPattern>> ret {};
	
	for(size_t i = 0, size = left.size(); i < size; ++i)
	{
		auto arg = dynamic_cast<Argument const*>(left[i].get());
		if (arg) {
			ret.first = i;
			ret.second = std::make_shared<Argument>(name(), arg->getValue());
			break;
		}
	}
	
	return ret;
}

std::pair<size_t, std::shared_ptr<LeafPattern>> Command::single_match(PatternList const& left) const
{
	std::pair<size_t, std::shared_ptr<LeafPattern>> ret {};
	
	for(size_t i = 0, size = left.size(); i < size; ++i)
	{
		auto arg = dynamic_cast<Argument const*>(left[i].get());
		if (arg) {
			if (name() == arg->getValue()) {
				ret.first = i;
				ret.second = std::make_shared<Command>(name(), value{true});
			}
			break;
		}
	}
	
	return ret;
}

std::pair<size_t, std::shared_ptr<LeafPattern>> Option::single_match(PatternList const& left) const
{
	std::pair<size_t, std::shared_ptr<LeafPattern>> ret {};
	
	for(size_t i = 0, size = left.size(); i < size; ++i)
	{
		auto leaf = std::dynamic_pointer_cast<LeafPattern>(left[i]);
		if (leaf && name() == leaf->name()) {
			ret.first = i;
			ret.second = leaf;
			break;
		}
	}
	
	return ret;
}

#pragma mark -
#pragma mark Parsing stuff

std::vector<PatternList> transform(PatternList pattern);

void BranchPattern::fix_repeating_arguments()
{
	std::vector<PatternList> either = transform(children());
	for(auto const& group : either) {
		// use multiset to help identify duplicate entries
		std::unordered_multiset<std::shared_ptr<Pattern>, PatternHasher> group_set {group.begin(), group.end()};
		for(auto const& e : group_set) {
			if (group_set.count(e) == 1)
				continue;
			
			LeafPattern* leaf = dynamic_cast<LeafPattern*>(e.get());
			if (!leaf) continue;
			
			bool ensureList = false;
			bool ensureInt = false;
			
			if (dynamic_cast<Command*>(leaf)) {
				ensureInt = true;
			} else if (dynamic_cast<Argument*>(leaf)) {
				ensureList = true;
			} else if (Option* o = dynamic_cast<Option*>(leaf)) {
				if (o->argCount()) {
					ensureList = true;
				} else {
					ensureInt = true;
				}
			}
			
			if (ensureList) {
				std::vector<std::string> newValue;
				if (leaf->getValue().isString()) {
					newValue = split(leaf->getValue().asString());
				}
				if (!leaf->getValue().isStringList()) {
					leaf->setValue(value{newValue});
				}
			} else if (ensureInt) {
				leaf->setValue(value{0});
			}
		}
	}
}

std::vector<PatternList> transform(PatternList pattern)
{
	std::vector<PatternList> result;
	
	std::vector<PatternList> groups;
	groups.emplace_back(std::move(pattern));
	
	while(!groups.empty()) {
		// pop off the first element
		auto children = std::move(groups[0]);
		groups.erase(groups.begin());
		
		// find the first branch node in the list
		auto child_iter = std::find_if(children.begin(), children.end(), [](std::shared_ptr<Pattern> const& p) {
			return dynamic_cast<BranchPattern const*>(p.get());
		});
		
		// no branch nodes left : expansion is complete for this grouping
		if (child_iter == children.end()) {
			result.emplace_back(std::move(children));
			continue;
		}
		
		// pop the child from the list
		auto child = std::move(*child_iter);
		children.erase(child_iter);
		
		// expand the branch in the appropriate way
		if (Either* either = dynamic_cast<Either*>(child.get())) {
			// "[e] + children" for each child 'e' in Either
			for(auto const& eitherChild : either->children()) {
				PatternList group = { eitherChild };
				group.insert(group.end(), children.begin(), children.end());
				
				groups.emplace_back(std::move(group));
			}
		} else if (OneOrMore* oneOrMore = dynamic_cast<OneOrMore*>(child.get())) {
			// child.children * 2 + children
			auto const& subchildren = oneOrMore->children();
			PatternList group = subchildren;
			group.insert(group.end(), subchildren.begin(), subchildren.end());
			group.insert(group.end(), children.begin(), children.end());
			
			groups.emplace_back(std::move(group));
		} else { // Required, Optional, OptionsShortcut
			BranchPattern* branch = dynamic_cast<BranchPattern*>(child.get());
			
			// child.children + children
			PatternList group = branch->children();
			group.insert(group.end(), children.begin(), children.end());
			
			groups.emplace_back(std::move(group));
		}
	}
	
	return result;
}

class Tokens {
public:
	Tokens(std::vector<std::string> tokens, bool isParsingArgv = true)
	: fTokens(std::move(tokens)),
	  fIsParsingArgv(isParsingArgv)
	{}
	
	explicit operator bool() const {
		return fIndex < fTokens.size();
	}

    static void acquire_if_nonempty(std::string *str, std::vector<std::string> *tokens) {
        if (! str->empty()) {
            tokens->resize(tokens->size() + 1);
            tokens->back().swap(*str);
            // str is now empty!
        }
    }
        
    static Tokens from_pattern(const std::string &source) {
        /* The following are tokens:
           [, ], (, ), |, ..., strings_without_spaces, <stuff in brackets>
           A string_without_spaces abutting <stuff_in_brackets> like so:
               foo<bar baz>
           count as only a single token.
        */
        
        bool within_brackets = false;
		std::vector<std::string> tokens;
        
        // When working on a string token, we append to this
        std::string currentStringToken;
        
        const char *sc = source.c_str();
        const size_t size = source.size();

        for (size_t i=0; i < size;)
        {
            const char c = sc[i];
            if (! within_brackets)
            {
                if (strchr("()[]|", c)) {
                    // Some sort of parenthesis
                    acquire_if_nonempty(&currentStringToken, &tokens);
                    tokens.emplace_back(1, c);
                    i += 1;
                } else if (strchr(" \t\n", c)) {
                    // Whitespace
                    acquire_if_nonempty(&currentStringToken, &tokens);
                    i += 1;
                } else if (! strncmp(sc + i, "...", 3)) {
                    // Ellipsis
                    acquire_if_nonempty(&currentStringToken, &tokens);
                    tokens.emplace_back("...");
                    i += strlen("...");
                } else {
                    // Normal char, or possibly bracket
                    currentStringToken.push_back(c);
                    within_brackets = (c == '<');
                    i += 1;
                }
            } else {
                // Within brackets
                currentStringToken.push_back(c);
                within_brackets = (c != '>');
                i += 1;
            }
        }
        // Grab trailing token
        acquire_if_nonempty(&currentStringToken, &tokens);
        
        return Tokens(tokens, false);
    }
	
	std::string const& current() const {
		if (*this)
			return fTokens[fIndex];
		
		static std::string const empty;
		return empty;
	}
	
	std::string the_rest() const {
		if (!*this)
			return {};
		return join(fTokens.begin()+fIndex,
			    fTokens.end(),
			    " ");
	}
	
	std::string pop() {
		return std::move(fTokens.at(fIndex++));
	}
	
	bool isParsingArgv() const { return fIsParsingArgv; }
	
	struct OptionError : std::runtime_error { using runtime_error::runtime_error; };
	
private:
	std::vector<std::string> fTokens;
	size_t fIndex = 0;
	bool fIsParsingArgv;
};
		
// Get all instances of 'T' from the pattern
template <typename T>
std::vector<T*> flat_filter(Pattern& pattern) {
	std::vector<Pattern*> flattened = pattern.flat([](Pattern const* p) -> bool {
		return dynamic_cast<T const*>(p);
	});
	
	// now, we're guaranteed to have T*'s, so just use static_cast
	std::vector<T*> ret;
	std::transform(flattened.begin(), flattened.end(), std::back_inserter(ret), [](Pattern* p) {
		return static_cast<T*>(p);
	});
	return ret;
}
        
/* Helper function to efficiently iterate over lines of a string 'source'. On input, line_end should be initialized to the start point for the iteration. On return, line_start will point at the next line, and line_end will point just after the trailing newline, or possibly at source.size(). The length of the line is line_end - line_start (and is guaranteed to be positive). Returns true if a line was returned, false if we reached the end. */
static bool get_next_line(const std::string &source, size_t *line_start, size_t *line_end)
{
    bool success = false;
    if (*line_end < source.size())
    {
        // Start at the end of the last line, or zero if this is the first call
        *line_start = *line_end;
        size_t newline = source.find('\n', *line_start);
        if (newline == std::string::npos) {
            // Point just after the last character
            *line_end = source.size();
        } else {
            // Point just after the newline
            *line_end = newline + 1;
        }
        // Empty lines are impossible
        assert(*line_end > *line_start);
        success = true;
    }
    return success;
}

/* Given a name like 'options:' and the usage text (as source), return a list of substrings of the
  source text that match the usage. The Python docopt is quite permissive, and is prone to false
  positive. For example, if the word "options:" appears anywhere in a line, it may be picked up. We
  wish to be more conservative. We model the source text as a list of lines, separated into sections.
  A section has a header like 'foo options:' that must not be indented (no leading whitespace) and must
  have a colon. All lines after the header until the next header are part of that section.
*/
std::vector<std::string> parse_section(const char *name, std::string const& source) {
    std::vector<std::string> sections;
    
    bool in_desired_section = false;
    size_t line_start = 0, line_end = 0;
    while (get_next_line(source, &line_start, &line_end)) {
        // It's a header line if the first character is not whitespace
        if (! isspace(source.at(line_start))) {
            // Check to see if the name is found before the first colon
            // Note that if name is not found at all, name_pos will have value npos, which is huge (and therefore not smaller than line_end)
            size_t name_pos = find_case_insensitive(source, name, line_start);
            in_desired_section = (name_pos < line_end && name_pos < source.find(':', line_start));
            if (in_desired_section) {
                // This line is the start of a section we want.
                // Append a blank line to our result vector to hold this section.
                sections.push_back(std::string());
            }
        }
        
        if (in_desired_section) {
            // We're in the section we want. Append the line to the current section.
            // Note this line may be its header.
            // Also note that result must be nonempty, because we always append to result at the point that in_desired_section is set to true
            sections.back().append(source, line_start, line_end - line_start);
        }
        
        // Update our next line to start just past the newline (or past the end of the source)
        line_start = line_end;
    }
    
    return sections;
}

bool is_argument_spec(std::string const& token) {
	if (token.empty())
		return false;
	
	if (token[0]=='<' && token[token.size()-1]=='>')
		return true;
	
	if (std::all_of(token.begin(), token.end(), &::isupper))
		return true;
	
	return false;
}

template <typename I>
std::vector<std::string> longOptions(I iter, I end) {
	std::vector<std::string> ret;
	std::transform(iter, end,
		       std::back_inserter(ret),
		       [](decltype(*iter)& opt) { return opt->longOption(); });
	return ret;
}

PatternList parse_long(Tokens& tokens, std::vector<Option>& options)
{
	// long ::= '--' chars [ ( ' ' | '=' ) chars ] ;
	std::string longOpt, equal;
	value val;
	std::tie(longOpt, equal, val) = partition(tokens.pop(), "=");
	
	assert(starts_with(longOpt, "--"));
	
	if (equal.empty()) {
		val = value{};
	}
	
	// detect with options match this long option
	std::vector<Option const*> similar;
	for(auto const& option : options) {
		if (option.longOption()==longOpt)
			similar.push_back(&option);
	}
	
	// maybe allow similar options that match by prefix
	if (tokens.isParsingArgv() && similar.empty()) {
		for(auto const& option : options) {
			if (option.longOption().empty())
				continue;
			if (starts_with(option.longOption(), longOpt))
				similar.push_back(&option);
		}
	}
	
	PatternList ret;
	
	if (similar.size() > 1) { // might be simply specified ambiguously 2+ times?
		std::vector<std::string> prefixes = longOptions(similar.begin(), similar.end());
		std::string error = "'" + longOpt + "' is not a unique prefix: ";
		error.append(join(prefixes.begin(), prefixes.end(), ", "));
		throw Tokens::OptionError(std::move(error));
	} else if (similar.empty()) {
		int argcount = equal.empty() ? 0 : 1;
		options.emplace_back("", longOpt, argcount);
		
		auto o = std::make_shared<Option>(options.back());
		if (tokens.isParsingArgv()) {
			o->setValue(argcount ? value{val} : value{true});
		}
		ret.push_back(o);
	} else {
		auto o = std::make_shared<Option>(*similar[0]);
		if (o->argCount() == 0) {
			if (val) {
				std::string error = o->longOption() + " must not have an argument";
				throw Tokens::OptionError(std::move(error));
			}
		} else {
			if (!val) {
				auto const& token = tokens.current();
				if (token.empty() || token=="--") {
					std::string error = o->longOption() + " requires an argument";
					throw Tokens::OptionError(std::move(error));
				}
				val = tokens.pop();
			}
		}
		if (tokens.isParsingArgv()) {
			o->setValue(val ? std::move(val) : value{true});
		}
		ret.push_back(o);
	}
	
	return ret;
}

PatternList parse_short(Tokens& tokens, std::vector<Option>& options)
{
	// shorts ::= '-' ( chars )* [ [ ' ' ] chars ] ;
	
	auto token = tokens.pop();
	
	assert(starts_with(token, "-"));
	assert(!starts_with(token, "--"));
	
	auto i = token.begin();
	++i; // skip the leading '-'
	
	PatternList ret;
	while (i != token.end()) {
		std::string shortOpt = { '-', *i };
		++i;
		
		std::vector<Option const*> similar;
		for(auto const& option : options) {
			if (option.shortOption()==shortOpt)
				similar.push_back(&option);
		}
		
		if (similar.size() > 1) {
			std::string error = shortOpt + " is specified ambiguously "
			+ std::to_string(similar.size()) + " times";
			throw Tokens::OptionError(std::move(error));
		} else if (similar.empty()) {
			options.emplace_back(shortOpt, "", 0);
			
			auto o = std::make_shared<Option>(options.back());
			if (tokens.isParsingArgv()) {
				o->setValue(value{true});
			}
			ret.push_back(o);
		} else {
			auto o = std::make_shared<Option>(*similar[0]);
			value val;
			if (o->argCount()) {
				if (i == token.end()) {
					// consume the next token
					auto const& token = tokens.current();
					if (token.empty() || token=="--") {
						std::string error = shortOpt + " requires an argument";
						throw Tokens::OptionError(std::move(error));
					}
					val = tokens.pop();
				} else {
					// consume all the rest
					val = std::string{i, token.end()};
					i = token.end();
				}
			}
			
			if (tokens.isParsingArgv()) {
				o->setValue(val ? std::move(val) : value{true});
			}
			ret.push_back(o);
		}
	}
	
	return ret;
}

PatternList parse_expr(Tokens& tokens, std::vector<Option>& options);

PatternList parse_atom(Tokens& tokens, std::vector<Option>& options)
{
	// atom ::= '(' expr ')' | '[' expr ']' | 'options'
	//             | long | shorts | argument | command ;
	
	std::string const& token = tokens.current();
	
	PatternList ret;
	
	if (token == "[") {
		tokens.pop();
		
		auto expr = parse_expr(tokens, options);
		
		auto trailing = tokens.pop();
		if (trailing != "]") {
			throw DocoptLanguageError("Mismatched '['");
		}
		
		ret.emplace_back(std::make_shared<Optional>(std::move(expr)));
	} else if (token=="(") {
		tokens.pop();
		
		auto expr = parse_expr(tokens, options);
		
		auto trailing = tokens.pop();
		if (trailing != ")") {
			throw DocoptLanguageError("Mismatched '('");
		}
		
		ret.emplace_back(std::make_shared<Required>(std::move(expr)));
	} else if (token == "options") {
		tokens.pop();
		ret.emplace_back(std::make_shared<OptionsShortcut>());
	} else if (starts_with(token, "--") && token != "--") {
		ret = parse_long(tokens, options);
	} else if (starts_with(token, "-") && token != "-" && token != "--") {
		ret = parse_short(tokens, options);
	} else if (is_argument_spec(token)) {
		ret.emplace_back(std::make_shared<Argument>(tokens.pop()));
	} else {
		ret.emplace_back(std::make_shared<Command>(tokens.pop()));
	}
	
	return ret;
}

PatternList parse_seq(Tokens& tokens, std::vector<Option>& options)
{
	// seq ::= ( atom [ '...' ] )* ;"""
	
	PatternList ret;
	
	while (tokens) {
		auto const& token = tokens.current();
		
		if (token=="]" || token==")" || token=="|")
			break;
		
		auto atom = parse_atom(tokens, options);
		if (tokens.current() == "...") {
			ret.emplace_back(std::make_shared<OneOrMore>(std::move(atom)));
			tokens.pop();
		} else {
			std::move(atom.begin(), atom.end(), std::back_inserter(ret));
		}
	}
	
	return ret;
}

std::shared_ptr<Pattern> maybe_collapse_to_required(PatternList&& seq)
{
	if (seq.size()==1) {
		return std::move(seq[0]);
	}
	return std::make_shared<Required>(std::move(seq));
}

std::shared_ptr<Pattern> maybe_collapse_to_either(PatternList&& seq)
{
	if (seq.size()==1) {
		return std::move(seq[0]);
	}
	return std::make_shared<Either>(std::move(seq));
}

PatternList parse_expr(Tokens& tokens, std::vector<Option>& options)
{
	// expr ::= seq ( '|' seq )* ;
	
	auto seq = parse_seq(tokens, options);
	
	if (tokens.current() != "|")
		return seq;
	
	PatternList ret;
	ret.emplace_back(maybe_collapse_to_required(std::move(seq)));
	
	while (tokens.current() == "|") {
		tokens.pop();
		seq = parse_seq(tokens, options);
		ret.emplace_back(maybe_collapse_to_required(std::move(seq)));
	}
	
	return { maybe_collapse_to_either(std::move(ret)) };
}

Required parse_pattern(std::string const& source, std::vector<Option>& options)
{
	auto tokens = Tokens::from_pattern(source);
	auto result = parse_expr(tokens, options);
	
	if (tokens)
		throw DocoptLanguageError("Unexpected ending: '" + tokens.the_rest() + "'");
	
	assert(result.size() == 1  &&  "top level is always one big");
	return Required{ std::move(result) };
}


std::string formal_usage(std::string const& section) {
	std::string ret = "(";
	
	auto i = section.find(':')+1;  // skip past "usage:"
	auto parts = split(section, i);
	for(size_t i = 1; i < parts.size(); ++i) {
		if (parts[i] == parts[0]) {
			ret += " ) | (";
		} else {
			ret.push_back(' ');
			ret += parts[i];
		}
	}
	
	ret += " )";
	return ret;
}

PatternList parse_argv(Tokens tokens, std::vector<Option>& options, bool options_first)
{
	// Parse command-line argument vector.
	//
	// If options_first:
	//    argv ::= [ long | shorts ]* [ argument ]* [ '--' [ argument ]* ] ;
	// else:
	//    argv ::= [ long | shorts | argument ]* [ '--' [ argument ]* ] ;
	
	PatternList ret;
	while (tokens) {
		auto const& token = tokens.current();
		
		if (token=="--") {
			// option list is done; convert all the rest to arguments
			while (tokens) {
				ret.emplace_back(std::make_shared<Argument>("", tokens.pop()));
			}
		} else if (starts_with(token, "--")) {
			auto&& parsed = parse_long(tokens, options);
			std::move(parsed.begin(), parsed.end(), std::back_inserter(ret));
		} else if (token[0]=='-' && token != "-") {
			auto&& parsed = parse_short(tokens, options);
			std::move(parsed.begin(), parsed.end(), std::back_inserter(ret));
		} else if (options_first) {
			// option list is done; convert all the rest to arguments
			while (tokens) {
				ret.emplace_back(std::make_shared<Argument>("", tokens.pop()));
			}
		} else {
			ret.emplace_back(std::make_shared<Argument>("", tokens.pop()));
		}
	}
	
	return ret;
}
    
std::vector<Option> parse_defaults(std::string const& doc) {
    /* Given the options section(s), extract the Option values from them
       An option may span lines like this:
     
           -v, --verbose: Be more verbose
                          and stuff like that
     
       We grab lines that have leading whitespace followed by a hyphen, up to
       the next such line.
    */
    std::vector<Option> defaults;
    
    for (const std::string &option_src : parse_section("options:", doc)) {
        // We walk over the defaults by line, building up an option into current_opt
        std::string current_opt;

         // Skip past "options:" by initializing line_end to just past the colon
        size_t option_label_pos = option_src.find(':');
        assert(option_label_pos != std::string::npos);
        size_t line_start = 0, line_end = option_label_pos + 1;
        
        // Walk over what's left by lines
        while (get_next_line(option_src, &line_start, &line_end)) {
            // See if this line has leading whitespace followed by a dash
            // The first non-whitespace character must be past line_start (so that we have leading whitespace),
            // and before line_end (so it's on this line)
            size_t first_non_whitespace = option_src.find_first_not_of(" \t", line_start);
            if (first_non_whitespace > line_start && first_non_whitespace < line_end && option_src.at(first_non_whitespace) == '-') {
                // It's a new option. Grab the current one.
                if (! current_opt.empty()) {
                    defaults.emplace_back(Option::parse(current_opt));
                }
                // Skip the leading whitespace in the option by starting at first_non_whitespace
                current_opt.assign(option_src, first_non_whitespace, line_end - first_non_whitespace);
            } else {
                // It's part of the option from the previous line. Append to it.
                if (current_opt.empty()) {
                    // TODO: Emit an error. This indicates that there was no leading dash on the first line, for example:
                    //    options: stuff
                } else {
                    current_opt.append(option_src, line_start, line_end);
                }
            }
        }
        // We've exhausted the lines; append the last option.
        if (! current_opt.empty()) {
            defaults.emplace_back(Option::parse(current_opt));
        }
    }
    return defaults;
}

bool isOptionSet(PatternList const& options, std::string const& opt1, std::string const& opt2 = "") {
	return std::any_of(options.begin(), options.end(), [&](std::shared_ptr<Pattern const> const& opt) -> bool {
		auto const& name = opt->name();
		if (name==opt1 || (!opt2.empty() && name==opt2)) {
			return opt->hasValue();
		}
		return false;
	});
}

void extras(bool help, bool version, PatternList const& options) {
	if (help && isOptionSet(options, "-h", "--help")) {
		throw DocoptExitHelp();
	}
	
	if (version && isOptionSet(options, "--version")) {
		throw DocoptExitVersion();
	}
}

// Parse the doc string and generate the Pattern tree
std::pair<Required, std::vector<Option>> create_pattern_tree(std::string const& doc)
{
	auto usage_sections = parse_section("usage:", doc);
	if (usage_sections.empty()) {
		throw DocoptLanguageError("'usage:' (case-insensitive) not found.");
	}
	if (usage_sections.size() > 1) {
		throw DocoptLanguageError("More than one 'usage:' (case-insensitive).");
	}
	
	std::vector<Option> options = parse_defaults(doc);
	Required pattern = parse_pattern(formal_usage(usage_sections[0]), options);
	
	std::vector<Option const*> pattern_options = flat_filter<Option const>(pattern);
	
	using UniqueOptions = std::unordered_set<Option const*, PatternHasher, PatternPointerEquality>;
	UniqueOptions const uniq_pattern_options { pattern_options.begin(), pattern_options.end() };
	
	// Fix up any "[options]" shortcuts with the actual option tree
	for(auto& options_shortcut : flat_filter<OptionsShortcut>(pattern)) {
		std::vector<Option> doc_options = parse_defaults(doc);
		
		// set(doc_options) - set(pattern_options)
		UniqueOptions uniq_doc_options;
		for(auto const& opt : doc_options) {
			if (uniq_pattern_options.count(&opt))
				continue;
			uniq_doc_options.insert(&opt);
		}
		
		// turn into shared_ptr's and set as children
		PatternList children;
		std::transform(uniq_doc_options.begin(), uniq_doc_options.end(),
			       std::back_inserter(children), [](Option const* opt) {
				       return std::make_shared<Option>(*opt);
			       });
		options_shortcut->setChildren(std::move(children));
	}
	
	return { std::move(pattern), std::move(options) };
}

std::map<std::string, value>
docopt::docopt_parse(std::string const& doc,
		     std::vector<std::string> const& argv,
		     bool help,
		     bool version,
		     bool options_first)
{
    bool log = false;
	Required pattern;
	std::vector<Option> options;
	try {
		std::tie(pattern, options) = create_pattern_tree(doc);
	} catch (Tokens::OptionError const& error) {
		throw DocoptLanguageError(error.what());
	}
	
	PatternList argv_patterns;
	try {
		argv_patterns = parse_argv(Tokens(argv), options, options_first);
	} catch (Tokens::OptionError const& error) {
		throw DocoptArgumentError(error.what());
	}
	
    if (log) {
        std::cerr << "Pattern: " << dump_pattern(pattern) << std::endl;
    }
    
	extras(help, version, argv_patterns);
	
	std::vector<std::shared_ptr<LeafPattern>> collected;
	bool matched = pattern.fix().match(argv_patterns, collected);
	if (matched && argv_patterns.empty()) {
		std::map<std::string, value> ret;
		
		// (a.name, a.value) for a in (pattern.flat() + collected)
		for (auto* p : pattern.leaves()) {
			ret[p->name()] = p->getValue();
		}
		
		for (auto const& p : collected) {
			ret[p->name()] = p->getValue();
		}
		
		return ret;
	}
	
	if (matched) {
		std::string leftover = join(argv.begin(), argv.end(), ", ");
		throw DocoptArgumentError("Unexpected argument: " + leftover);
	}
	
	throw DocoptArgumentError("Arguments did not match expected patterns"); // BLEH. Bad error.
}
		
std::map<std::string, value>
docopt::docopt(std::string const& doc,
	       std::vector<std::string> const& argv,
	       bool help,
	       std::string const& version,
	       bool options_first) noexcept
{
	try {
		return docopt_parse(doc, argv, help, !version.empty(), options_first);
	} catch (DocoptExitHelp const&) {
		std::cout << doc << std::endl;
		std::exit(0);
	} catch (DocoptExitVersion const&) {
		std::cout << version << std::endl;
		std::exit(0);
	} catch (DocoptLanguageError const& error) {
		std::cerr << "Docopt usage string could not be parsed" << std::endl;
		std::cerr << error.what() << std::endl;
		std::exit(-1);
	} catch (DocoptArgumentError const& error) {
		std::cerr << error.what();
		std::cout << std::endl;
		std::cout << doc << std::endl;
		std::exit(-1);
	} /* Any other exception is unexpected: let std::terminate grab it */
}
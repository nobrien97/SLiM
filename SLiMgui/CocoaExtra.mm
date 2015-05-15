//
//  CocoaExtra.m
//  SLiM
//
//  Created by Ben Haller on 1/22/15.
//  Copyright (c) 2015 Philipp Messer.  All rights reserved.
//	A product of the Messer Lab, http://messerlab.org/software/
//

//	This file is part of SLiM.
//
//	SLiM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
//	SLiM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License along with SLiM.  If not, see <http://www.gnu.org/licenses/>.


#import "CocoaExtra.h"

#include "script.h"

#include <stdexcept>


@implementation SLiMSyntaxColoredTextView

// produce standard text attributes including our font (Menlo 11), tab stops (every three spaces), and font colors for syntax coloring
+ (NSDictionary *)consoleTextAttributesWithColor:(NSColor *)textColor
{
	static NSFont *menlo11Font = nil;
	static NSMutableParagraphStyle *paragraphStyle = nil;
	
	if (!menlo11Font)
		menlo11Font = [[NSFont fontWithName:@"Menlo" size:11.0] retain];
	
	if (!paragraphStyle)
	{
		paragraphStyle = [[NSParagraphStyle defaultParagraphStyle] mutableCopy];
		
		CGFloat tabInterval = [menlo11Font maximumAdvancement].width * 3;
		NSMutableArray *tabs = [NSMutableArray array];
		
		[paragraphStyle setDefaultTabInterval:tabInterval];
		
		for (int tabStop = 1; tabStop <= 20; ++tabStop)
			[tabs addObject:[[NSTextTab alloc] initWithTextAlignment:NSLeftTextAlignment location:tabInterval * tabStop options:nil]];
		
		[paragraphStyle setTabStops:tabs];
	}
	
	if (textColor)
		return [NSDictionary dictionaryWithObjectsAndKeys:textColor, NSForegroundColorAttributeName, menlo11Font, NSFontAttributeName, paragraphStyle, NSParagraphStyleAttributeName, nil];
	else
		return [NSDictionary dictionaryWithObjectsAndKeys:menlo11Font, NSFontAttributeName, paragraphStyle, NSParagraphStyleAttributeName, nil];
}

// handle autoindent by matching the whitespace beginning the current line
- (void)insertNewline:(id)sender
{
	NSString *textString = [self string];
	NSUInteger selectionStart = [self selectedRange].location;
	NSCharacterSet *newlineChars = [NSCharacterSet newlineCharacterSet];
	NSCharacterSet *whitespaceChars = [NSCharacterSet whitespaceCharacterSet];
	
	// start at the start of the selection and move backwards to the beginning of the line
	NSUInteger lineStart = selectionStart;
	
	while (lineStart > 0)
	{
		unichar ch = [textString characterAtIndex:lineStart - 1];
		
		if ([newlineChars characterIsMember:ch])
			break;
		
		--lineStart;
	}
	
	// now we're either at the beginning of the content, or the beginning of the line; now find the end of the whitespace there, up to where we started
	NSUInteger whitespaceEnd = lineStart;
	
	while (whitespaceEnd < selectionStart)
	{
		unichar ch = [textString characterAtIndex:whitespaceEnd];
		
		if (![whitespaceChars characterIsMember:ch])
			break;
		
		++whitespaceEnd;
	}
	
	// now we have the range of the leading whitespace; copy that, call super to insert the newline, and then paste in the whitespace
	NSRange whitespaceRange = NSMakeRange(lineStart, whitespaceEnd - lineStart);
	NSString *whitespaceString = [textString substringWithRange:whitespaceRange];
	
	[super insertNewline:sender];
	[self insertText:whitespaceString];
}

// NSTextView copies only plain text for us, because it is set to have rich text turned off.  That setting only means it is turned off for the user; the
// user can't change the font, size, etc.  But we still can, and do, programatically to do our syntax formatting.  We want that style information to get
// copied to the pasteboard, and as far as I can tell this subclass is necessary to make it happen.  Seems kind of lame.
- (IBAction)copy:(id)sender
{
	NSPasteboard *pasteboard = [NSPasteboard generalPasteboard];
	NSAttributedString *attrString = [self textStorage];
	NSRange selectedRange = [self selectedRange];
	NSAttributedString *attrStringInRange = [attrString attributedSubstringFromRange:selectedRange];
	
	// The documentation sucks, but as far as I can tell, this puts both a plain-text and a rich-text representation on the pasteboard
	[pasteboard clearContents];
	[pasteboard writeObjects:[NSArray arrayWithObject:attrStringInRange]];
}

- (IBAction)shiftSelectionLeft:(id)sender
{
	if ([self isEditable])
	{
		NSTextStorage *ts = [self textStorage];
		NSMutableString *scriptString = [[self string] mutableCopy];
		int scriptLength = (int)[scriptString length];
		NSRange selectedRange = [self selectedRange];
		NSCharacterSet *newlineChars = [NSCharacterSet newlineCharacterSet];
		NSUInteger scanPosition;
		
		// start at the start of the selection and scan backwards over non-newline text until we hit a newline or the start of the file
		scanPosition = selectedRange.location;
		
		while (scanPosition > 0)
		{
			if ([newlineChars characterIsMember:[scriptString characterAtIndex:scanPosition - 1]])
				break;
			
			--scanPosition;
		}
		
		// ok, we're at the start of the line that the selection starts on; start removing tabs
		[ts beginEditing];
		
		while ((scanPosition == selectedRange.location) || (scanPosition < selectedRange.location + selectedRange.length))
		{
			// if we are at the very end of the script string, then we have hit the end and we're done
			if (scanPosition == scriptLength)
				break;
			
			// insert a tab at the start of this line and adjust our selection
			if ([scriptString characterAtIndex:scanPosition] == '\t')
			{
				[ts replaceCharactersInRange:NSMakeRange(scanPosition, 1) withString:@""];
				[scriptString replaceCharactersInRange:NSMakeRange(scanPosition, 1) withString:@""];
				scriptLength--;
				
				if (scanPosition < selectedRange.location)
					selectedRange.location--;
				else
					selectedRange.length--;
			}
			
			// now scan forward to the end of this line
			while (scanPosition < scriptLength)
			{
				if ([newlineChars characterIsMember:[scriptString characterAtIndex:scanPosition]])
					break;
				
				++scanPosition;
			}
			
			// and then scan forward to the beginning of the next line
			while (scanPosition < scriptLength)
			{
				if (![newlineChars characterIsMember:[scriptString characterAtIndex:scanPosition]])
					break;
				
				++scanPosition;
			}
		}
		
		[ts endEditing];
		[self setSelectedRange:selectedRange];
	}
	else
	{
		NSBeep();
	}
}

- (IBAction)shiftSelectionRight:(id)sender
{
	if ([self isEditable])
	{
		NSTextStorage *ts = [self textStorage];
		NSMutableString *scriptString = [[self string] mutableCopy];
		int scriptLength = (int)[scriptString length];
		NSRange selectedRange = [self selectedRange];
		NSCharacterSet *newlineChars = [NSCharacterSet newlineCharacterSet];
		NSUInteger scanPosition;
		
		// start at the start of the selection and scan backwards over non-newline text until we hit a newline or the start of the file
		scanPosition = selectedRange.location;
		
		while (scanPosition > 0)
		{
			if ([newlineChars characterIsMember:[scriptString characterAtIndex:scanPosition - 1]])
				break;
			
			--scanPosition;
		}
		
		// ok, we're at the start of the line that the selection starts on; start inserting tabs
		[ts beginEditing];
		
		while ((scanPosition == selectedRange.location) || (scanPosition < selectedRange.location + selectedRange.length))
		{
			// insert a tab at the start of this line and adjust our selection
			[ts replaceCharactersInRange:NSMakeRange(scanPosition, 0) withString:@"\t"];
			[scriptString replaceCharactersInRange:NSMakeRange(scanPosition, 0) withString:@"\t"];
			scriptLength++;
			
			if ((scanPosition < selectedRange.location) || (selectedRange.length == 0))
				selectedRange.location++;
			else
				selectedRange.length++;
			
			// now scan forward to the end of this line
			while (scanPosition < scriptLength)
			{
				if ([newlineChars characterIsMember:[scriptString characterAtIndex:scanPosition]])
					break;
				
				++scanPosition;
			}
			
			// and then scan forward to the beginning of the next line
			while (scanPosition < scriptLength)
			{
				if (![newlineChars characterIsMember:[scriptString characterAtIndex:scanPosition]])
					break;
				
				++scanPosition;
			}
			
			// if we are at the very end of the script string, then we have hit the end and we're done
			if (scanPosition == scriptLength)
				break;
		}
		
		[ts endEditing];
		[self setSelectedRange:selectedRange];
	}
	else
	{
		NSBeep();
	}
}

// Check whether a token string is a special identifier like "pX", "gX", or "mX"
- (BOOL)tokenStringIsSpecialIdentifier:(const std::string &)token_string
{
	int len = (int)token_string.length();
	
	if (len >= 2)
	{
		unichar first_ch = token_string[0];
		
		if ((first_ch == 'p') || (first_ch == 'g') || (first_ch == 'm'))
		{
			for (int ch_index = 1; ch_index < len; ++ch_index)
			{
				unichar idx_ch = token_string[ch_index];
				
				if ((idx_ch < '0') || (idx_ch > '9'))
					return NO;
			}
			
			return YES;
		}
	}
	
	return NO;
}

- (void)syntaxColorForSLiMScript
{
	// Construct a Script object from the current script string
	NSString *scriptString = [self string];
	std::string script_string([scriptString UTF8String]);
	Script script(1, 1, script_string, 0);
	
	// Tokenize
	try
	{
		script.Tokenize(true);	// keep nonsignificant tokens - whitespace and comments
	}
	catch (std::runtime_error err)
	{
		// if we get a raise, we just use as many tokens as we got
		//NSString *errorString = [NSString stringWithUTF8String:GetUntrimmedRaiseMessage().c_str()];
		//NSLog(@"raise during syntax coloring tokenization: %@", errorString);
	}
	
	// Set up our shared colors
	static NSColor *numberLiteralColor = nil;
	static NSColor *stringLiteralColor = nil;
	static NSColor *commentColor = nil;
	static NSColor *identifierColor = nil;
	static NSColor *keywordColor = nil;
	
	if (!numberLiteralColor)
	{
		numberLiteralColor = [[NSColor colorWithCalibratedRed:28/255.0 green:0/255.0 blue:207/255.0 alpha:1.0] retain];
		stringLiteralColor = [[NSColor colorWithCalibratedRed:196/255.0 green:26/255.0 blue:22/255.0 alpha:1.0] retain];
		commentColor = [[NSColor colorWithCalibratedRed:0/255.0 green:116/255.0 blue:0/255.0 alpha:1.0] retain];
		identifierColor = [[NSColor colorWithCalibratedRed:63/255.0 green:110/255.0 blue:116/255.0 alpha:1.0] retain];
		keywordColor = [[NSColor colorWithCalibratedRed:170/255.0 green:13/255.0 blue:145/255.0 alpha:1.0] retain];
	}
	
	// Syntax color!
	NSTextStorage *ts = [self textStorage];
	
	[ts beginEditing];
	
	[ts removeAttribute:NSForegroundColorAttributeName range:NSMakeRange(0, [ts length])];
	
	for (ScriptToken *token : script.Tokens())
	{
		NSRange tokenRange = NSMakeRange(token->token_start_, token->token_end_ - token->token_start_ + 1);
		
		if (token->token_type_ == TokenType::kTokenNumber)
			[ts addAttribute:NSForegroundColorAttributeName value:numberLiteralColor range:tokenRange];
		if (token->token_type_ == TokenType::kTokenString)
			[ts addAttribute:NSForegroundColorAttributeName value:stringLiteralColor range:tokenRange];
		if (token->token_type_ == TokenType::kTokenComment)
			[ts addAttribute:NSForegroundColorAttributeName value:commentColor range:tokenRange];
		if (token->token_type_ > TokenType::kFirstIdentifierLikeToken)
			[ts addAttribute:NSForegroundColorAttributeName value:keywordColor range:tokenRange];
		if (token->token_type_ == TokenType::kTokenIdentifier)
		{
			// most identifiers are left as black; only special ones get colored
			const std::string &token_string = token->token_string_;
			
			if ((token_string.compare("T") == 0) ||
				(token_string.compare("F") == 0) ||
				(token_string.compare("E") == 0) ||
				(token_string.compare("PI") == 0) ||
				(token_string.compare("INF") == 0) ||
				(token_string.compare("NAN") == 0) ||
				(token_string.compare("NULL") == 0) ||
				(token_string.compare("sim") == 0) ||
				[self tokenStringIsSpecialIdentifier:token_string])
				[ts addAttribute:NSForegroundColorAttributeName value:identifierColor range:tokenRange];
		}
	}
	
	[ts endEditing];
}

- (void)syntaxColorForSLiMInput
{
	NSTextStorage *textStorage = [self textStorage];
	NSString *string = [self string];
	NSArray *lines = [string componentsSeparatedByString:@"\n"];
	int lineCount = (int)[lines count];
	int stringPosition = 0;
	
	// Set up our shared attributes
	static NSDictionary *poundDirectiveAttrs = nil;
	static NSDictionary *commentAttrs = nil;
	static NSDictionary *subpopAttrs = nil;
	static NSDictionary *genomicElementAttrs = nil;
	static NSDictionary *mutationTypeAttrs = nil;
	
	if (!poundDirectiveAttrs)
	{
		poundDirectiveAttrs = [[NSDictionary alloc] initWithObjectsAndKeys:[NSColor colorWithCalibratedRed:196/255.0 green:26/255.0 blue:22/255.0 alpha:1.0], NSForegroundColorAttributeName, nil];
		commentAttrs = [[NSDictionary alloc] initWithObjectsAndKeys:[NSColor colorWithCalibratedRed:0/255.0 green:116/255.0 blue:0/255.0 alpha:1.0], NSForegroundColorAttributeName, nil];
		subpopAttrs = [[NSDictionary alloc] initWithObjectsAndKeys:[NSColor colorWithCalibratedRed:28/255.0 green:0/255.0 blue:207/255.0 alpha:1.0], NSForegroundColorAttributeName, nil];
		genomicElementAttrs = [[NSDictionary alloc] initWithObjectsAndKeys:[NSColor colorWithCalibratedRed:63/255.0 green:110/255.0 blue:116/255.0 alpha:1.0], NSForegroundColorAttributeName, nil];
		mutationTypeAttrs = [[NSDictionary alloc] initWithObjectsAndKeys:[NSColor colorWithCalibratedRed:170/255.0 green:13/255.0 blue:145/255.0 alpha:1.0], NSForegroundColorAttributeName, nil];
	}
	
	// And then tokenize and color
	[textStorage beginEditing];
	[textStorage removeAttribute:NSForegroundColorAttributeName range:NSMakeRange(0, [textStorage length])];
	
	for (int lineIndex = 0; lineIndex < lineCount; ++lineIndex)
	{
		NSString *line = [lines objectAtIndex:lineIndex];
		NSRange lineRange = NSMakeRange(stringPosition, (int)[line length]);
		int nextStringPosition = (int)(stringPosition + lineRange.length + 1);			// +1 for the newline
		
		if (lineRange.length)
		{
			//NSLog(@"lineIndex %d, lineRange == %@", lineIndex, NSStringFromRange(lineRange));
			
			// find comments and color and remove them
			NSRange commentRange = [line rangeOfString:@"//"];
			
			if ((commentRange.location != NSNotFound) && (commentRange.length == 2))
			{
				int commentLength = (int)(lineRange.length - commentRange.location);
				
				[textStorage addAttributes:commentAttrs range:NSMakeRange(lineRange.location + commentRange.location, commentLength)];
				
				lineRange.length -= commentLength;
				line = [line substringToIndex:commentRange.location];
			}
			
			// if anything is left...
			if (lineRange.length)
			{
				// remove leading whitespace
				do {
					NSRange leadingWhitespaceRange = [line rangeOfCharacterFromSet:[NSCharacterSet whitespaceCharacterSet] options:NSAnchoredSearch];
					
					if (leadingWhitespaceRange.location == NSNotFound || leadingWhitespaceRange.length == 0)
						break;
					
					lineRange.location += leadingWhitespaceRange.length;
					lineRange.length -= leadingWhitespaceRange.length;
					line = [line substringFromIndex:leadingWhitespaceRange.length];
				} while (YES);
				
				// remove trailing whitespace
				do {
					NSRange trailingWhitespaceRange = [line rangeOfCharacterFromSet:[NSCharacterSet whitespaceCharacterSet] options:NSAnchoredSearch | NSBackwardsSearch];
					
					if (trailingWhitespaceRange.location == NSNotFound || trailingWhitespaceRange.length == 0)
						break;
					
					lineRange.length -= trailingWhitespaceRange.length;
					line = [line substringToIndex:trailingWhitespaceRange.location];
				} while (YES);
				
				// if anything is left...
				if (lineRange.length)
				{
					// find pound directives and color them
					if ([line characterAtIndex:0] == '#')
						[textStorage addAttributes:poundDirectiveAttrs range:lineRange];
					else
					{
						NSRange scanRange = NSMakeRange(0, lineRange.length);
						
						do {
							NSRange tokenRange = [line rangeOfString:@"\\b[pgm][0-9]+\\b" options:NSRegularExpressionSearch range:scanRange];
							
							if (tokenRange.location == NSNotFound || tokenRange.length == 0)
								break;
							
							NSString *substring = [line substringWithRange:tokenRange];
							NSDictionary *syntaxAttrs = nil;
							
							if ([substring characterAtIndex:0] == 'p')
								syntaxAttrs = subpopAttrs;
							else if ([substring characterAtIndex:0] == 'g')
								syntaxAttrs = genomicElementAttrs;
							else if ([substring characterAtIndex:0] == 'm')
								syntaxAttrs = mutationTypeAttrs;
							
							if (syntaxAttrs)
								[textStorage addAttributes:syntaxAttrs range:NSMakeRange(tokenRange.location + lineRange.location, tokenRange.length)];
							
							scanRange.length = (scanRange.location + scanRange.length) - (tokenRange.location + tokenRange.length);
							scanRange.location = (tokenRange.location + tokenRange.length);
							
							if (scanRange.length < 2)
								break;
						} while (YES);
					}
				}
			}
		}
		
		stringPosition = nextStringPosition;
	}
	
	[textStorage endEditing];
}

- (void)clearSyntaxColoring
{
	NSTextStorage *textStorage = [self textStorage];
	
	[textStorage beginEditing];
	[textStorage removeAttribute:NSForegroundColorAttributeName range:NSMakeRange(0, [textStorage length])];
	[textStorage endEditing];
}

// add delegation of the -rangeForUserCompletion method so SLiMscribe can handle completion all by itself

- (NSRange)rangeForUserCompletion
{
	NSRange range = [super rangeForUserCompletion];
	
	if ([[self delegate] respondsToSelector:@selector(textView:rangeForUserCompletion:)])
		range = [(id <SLiMSyntaxColoredTextViewDelegate>)[self delegate] textView:self rangeForUserCompletion:range];
	
	return range;
}

@end


//
//	The rest of this file is live only in SLiMgui, not in SLiMscribe
//
#ifdef SLIMGUI

#import "AppDelegate.h"


@implementation SLiMTableView

- (BOOL)acceptsFirstResponder
{
	return NO;
}

@end


@implementation SLiMApp

- (SLiMWindowController *)mainSLiMWindowController
{
	return (SLiMWindowController *)[[self mainWindow] windowController];
}

- (IBAction)showHelp:(id)sender
{
	// Delegate this to, hey, our delegate.  Annoying that Apple hasn't fixed this so subclassing NSApplication is not necessary...
	[(AppDelegate *)[self delegate] showHelp:sender];
}

@end


@implementation SLiMDocumentController

@end


static NSDictionary *tickAttrs = nil;
static NSDictionary *disabledTickAttrs = nil;
const int numberOfTicks = 4;
const int tickLength = 5;
const int heightForTicks = 16;

@implementation SLiMColorStripeView

+ (void)initialize
{
	if (!tickAttrs)
		tickAttrs = [[NSDictionary alloc] initWithObjectsAndKeys:[NSColor blackColor], NSForegroundColorAttributeName, [NSFont systemFontOfSize:9.0], NSFontAttributeName, nil];
	if (!disabledTickAttrs)
		disabledTickAttrs = [[NSDictionary alloc] initWithObjectsAndKeys:[NSColor colorWithCalibratedWhite:0.6 alpha:1.0], NSForegroundColorAttributeName, [NSFont systemFontOfSize:9.0], NSFontAttributeName, nil];
	
	[self exposeBinding:@"enabled"];
}

- (void)setEnabled:(BOOL)enabled
{
	if (_enabled != enabled)
	{
		_enabled = enabled;
		
		[self setNeedsDisplay:YES];
	}
}

- (void)awakeFromNib
{
	[self bind:@"enabled" toObject:[[self window] windowController] withKeyPath:@"invalidSimulation" options:[NSDictionary dictionaryWithObject:NSNegateBooleanTransformerName forKey:NSValueTransformerNameBindingOption]];
}

- (void)dealloc
{
	[self unbind:@"enabled"];
	
	[super dealloc];
}

- (BOOL)isOpaque
{
	return NO;	// we have a 10-pixel margin on our left and right to allow tick labels to overflow our apparent bounds
}

- (void)setMetricToPlot:(int)metricToPlot
{
	if (_metricToPlot != metricToPlot)
	{
		_metricToPlot = metricToPlot;
		
		[self setNeedsDisplay:YES];
	}
}

- (void)drawTicksInContentRect:(NSRect)contentRect
{
	BOOL enabled = [self enabled];
	int metric = [self metricToPlot];
	NSRect interiorRect = NSInsetRect(contentRect, 1, 1);
	
	for (int position = 0; position <= numberOfTicks; ++position)
	{
		double fraction = position / (double)numberOfTicks;
		int tickLeft = (int)floor(interiorRect.origin.x + fraction * (interiorRect.size.width - 1));
		NSRect tickRect = NSMakeRect(tickLeft, contentRect.origin.y - tickLength, 1, tickLength);
		
		[[NSColor colorWithCalibratedWhite:(enabled ? 0.5 : 0.6) alpha:1.0] set];
		NSRectFill(tickRect);
		
		NSString *tickLabel = nil;
		
		if (metric == 1)
		{
			if (position == 0) tickLabel = @"0.0";
			if (position == 1) tickLabel = @"0.5";
			if (position == 2) tickLabel = @"1.0";
			if (position == 3) tickLabel = @"2.0";
			if (position == 4) tickLabel = @"∞";
		}
		else if (metric == 2)
		{
			if (position == 0) tickLabel = @"−0.5";
			if (position == 1) tickLabel = @"−0.25";
			if (position == 2) tickLabel = @"0.0";
			if (position == 3) tickLabel = @"0.5";
			if (position == 4) tickLabel = @"1.0";
		}
		
		if (tickLabel)
		{
			NSAttributedString *tickAttrLabel = [[NSAttributedString alloc] initWithString:tickLabel attributes:(enabled ? tickAttrs : disabledTickAttrs)];
			NSSize tickLabelSize = [tickAttrLabel size];
			int tickLabelX = tickLeft - (int)round(tickLabelSize.width / 2.0);
			
			[tickAttrLabel drawAtPoint:NSMakePoint(tickLabelX, contentRect.origin.y - (tickLength + 12))];
			
			[tickAttrLabel release];
		}
	}
}

- (void)drawRect:(NSRect)dirtyRect
{
	NSRect bounds = [self bounds];
	NSRect contentRect = NSMakeRect(bounds.origin.x + 10, bounds.origin.y + heightForTicks, bounds.size.width - 20, bounds.size.height - heightForTicks);
	NSRect interiorRect = NSInsetRect(contentRect, 1, 1);
	
	// frame the content area itself
	[[NSColor colorWithCalibratedWhite:0.6 alpha:1.0] set];
	NSFrameRect(contentRect);
	
	double scaling = [self scalingFactor];
	int metric = [self metricToPlot];
	
	// draw our stripe
	for (int x = 0; x < interiorRect.size.width; ++x)
	{
		NSRect stripe = NSMakeRect(interiorRect.origin.x + x, interiorRect.origin.y, 1, interiorRect.size.height);
		float red = 0.0, green = 0.0, blue = 0.0;
		double fraction = x / interiorRect.size.width;	// note this is never quite 1.0
		
		if (metric == 1)
		{
			double fitness;
			
			if (fraction < 0.5) fitness = fraction * 2.0;						// [0.0, 0.5] -> [0.0, 1.0]
			else if (fraction < 0.75) fitness = (fraction - 0.5) * 4.0 + 1.0;	// [0.5, 0.75] -> [1.0, 2.0]
			else fitness = 0.50 / (0.25 - (fraction - 0.75));					// [0.75, 1.0] -> [2.0, +Inf]
			
			RGBForFitness(fitness, &red, &green, &blue, scaling);
		}
		else if (metric == 2)
		{
			double selectionCoeff;
			
			if (fraction < 0.5) selectionCoeff = fraction - 0.5;				// [0.0, 0.5] -> [-0.5, 0]
			else selectionCoeff = (fraction - 0.5) * 2.0;						// [0.5, 1.0] -> [0.0, 1.0]
			
			RGBForSelectionCoeff(selectionCoeff, &red, &green, &blue, scaling);
		}
		
		[[NSColor colorWithCalibratedRed:red green:green blue:blue alpha:1.0] set];
		NSRectFill(stripe);
	}
	
	// draw ticks at bottom of content rect
	[self drawTicksInContentRect:contentRect];
	
	// if we are not enabled, wash over the interior with light gray
	if (![self enabled])
	{
		[[NSColor colorWithCalibratedWhite:0.9 alpha:0.8] set];
		NSRectFillUsingOperation(interiorRect, NSCompositeSourceOver);
	}
}

@end

const float greenBrightness = 0.9f;

void RGBForFitness(double value, float *colorRed, float *colorGreen, float *colorBlue, double scalingFactor)
{
	// apply the scaling factor
	value = (value - 1.0) * scalingFactor + 1.0;
	
	if (value <= 0.5)
	{
		// value <= 0.5 is a shade of red, going down to black
		*colorRed = (float)(value * 2.0);
		*colorGreen = 0.0;
		*colorBlue = 0.0;
	}
	else if (value >= 2.0)
	{
		// value >= 2.0 is a shade of green, going up to white
		*colorRed = (float)((value - 2.0) * greenBrightness / value);
		*colorGreen = greenBrightness;
		*colorBlue = (float)((value - 2.0) * greenBrightness / value);
	}
	else if (value <= 1.0)
	{
		// value <= 1.0 (but > 0.5) goes from red (unfit) to yellow (neutral)
		*colorRed = 1.0;
		*colorGreen = (float)((value - 0.5) * 2.0);
		*colorBlue = 0.0;
	}
	else	// 1.0 < value < 2.0
	{
		// value > 1.0 (but < 2.0) goes from yellow (neutral) to green (fit)
		*colorRed = (float)(2.0 - value);
		*colorGreen = (float)(greenBrightness + (1.0 - greenBrightness) * (2.0 - value));
		*colorBlue = 0.0;
	}
}

void RGBForSelectionCoeff(double value, float *colorRed, float *colorGreen, float *colorBlue, double scalingFactor)
{
	// apply a scaling factor; this could be user-adjustible since different models have different relevant fitness ranges
	value *= scalingFactor;
	
	// and add 1, just so we can re-use the same code as in RGBForFitness()
	value += 1.0;
	
	if (value <= 0.5)
	{
		// value <= 0.5 is a shade of red, going down to black
		*colorRed = (float)(value * 2.0);
		*colorGreen = 0.0;
		*colorBlue = 0.0;
	}
	else if (value >= 2.0)
	{
		// value >= 2.0 is a shade of green, going up to white
		*colorRed = (float)((value - 2.0) * greenBrightness / value);
		*colorGreen = greenBrightness;
		*colorBlue = (float)((value - 2.0) * greenBrightness / value);
	}
	else if (value <= 1.0)
	{
		// value <= 1.0 (but > 0.5) goes from red (unfit) to yellow (neutral)
		*colorRed = 1.0;
		*colorGreen = (float)((value - 0.5) * 2.0);
		*colorBlue = 0.0;
	}
	else	// 1.0 < value < 2.0
	{
		// value > 1.0 (but < 2.0) goes from yellow (neutral) to green (fit)
		*colorRed = (float)(2.0 - value);
		*colorGreen = (float)(greenBrightness + (1.0 - greenBrightness) * (2.0 - value));
		*colorBlue = 0.0;
	}
}


@implementation SLiMWhiteView

- (void)drawRect:(NSRect)dirtyRect
{
	[[NSColor whiteColor] set];
	NSRectFill(dirtyRect);
}

@end


@implementation SLiMMenuButton

- (void)fixMenu
{
	NSMenu *menu = [self slimMenu];
	NSDictionary *itemAttrs = [NSDictionary dictionaryWithObject:[NSFont systemFontOfSize:12.0] forKey:NSFontAttributeName];
	
	for (int i = 0; i < [menu numberOfItems]; ++i)
	{
		NSMenuItem *menuItem = [menu itemAtIndex:i];
		NSString *title = [menuItem title];
		NSAttributedString *attrTitle = [[NSAttributedString alloc] initWithString:title attributes:itemAttrs];
		
		[menuItem setAttributedTitle:attrTitle];
		[attrTitle release];
	}
}

- (void)mouseDown:(NSEvent *)theEvent
{
	// We do not call super; we do mouse tracking entirely ourselves
	//[super mouseDown:theEvent];
	
	if ([self isEnabled])
	{
		NSRect bounds = [self bounds];
		
		[self highlight:YES];
		
		[self fixMenu];
		[[self slimMenu] popUpMenuPositioningItem:nil atLocation:NSMakePoint(bounds.size.width * 0.80 - 1, bounds.size.height * 0.80 + 1) inView:self];
		
		[self highlight:NO];
	}
}

- (void)mouseDragged:(NSEvent *)theEvent
{
	// We do not call super; we do mouse tracking entirely ourselves
	//[super mouseDragged:theEvent];
}

- (void)mouseUp:(NSEvent *)theEvent
{
	// We do not call super; we do mouse tracking entirely ourselves
	//[super mouseDown:theEvent];
}

@end

@implementation SLiMColorCell

- (void)setObjectValue:(id)objectValue
{
	// Ensure that only NSColor objects are set as our object value
	if ([objectValue isKindOfClass:[NSColor class]])
		[super setObjectValue:objectValue];
}

- (void)drawWithFrame:(NSRect)cellFrame inView:(NSView *)controlView
{
	NSRect swatchFrame = NSInsetRect(cellFrame, 1, 1);
	
	// Make our swatch be square, if we have enough room to do it
	if (swatchFrame.size.width > swatchFrame.size.height)
	{
		swatchFrame.origin.x = (int)round(swatchFrame.origin.x + (swatchFrame.size.width - swatchFrame.size.height) / 2.0);
		swatchFrame.size.width = swatchFrame.size.height;
	}
	
	[[NSColor blackColor] set];
	NSFrameRect(swatchFrame);
	
	[(NSColor *)[self objectValue] set];
	NSRectFill(NSInsetRect(swatchFrame, 1, 1));
}

@end


@implementation SLiMSelectionView

- (void)drawRect:(NSRect)dirtyRect
{
	static NSDictionary *labelAttrs = nil;
	
	if (!labelAttrs)
		labelAttrs = [[NSDictionary alloc] initWithObjectsAndKeys:[NSFont fontWithName:@"Times New Roman" size:10], NSFontAttributeName, [NSColor blackColor], NSForegroundColorAttributeName, nil];
	
	NSRect bounds = [self bounds];
	SLiMSelectionMarker *marker = (SLiMSelectionMarker *)[self window];
	NSAttributedString *attrLabel = [[NSAttributedString alloc] initWithString:[marker label] attributes:labelAttrs];
	NSSize labelStringSize = [attrLabel size];
	NSSize labelSize = NSMakeSize(round(labelStringSize.width + 8.0), round(labelStringSize.height + 1.0));
	NSRect labelRect = NSMakeRect(bounds.origin.x + round(bounds.size.width / 2.0), bounds.origin.y + bounds.size.height - labelSize.height, labelSize.width, labelSize.height);
	
	if ([marker isLeftMarker])
		labelRect.origin.x -= round(labelSize.width - 1.0);
	
	// Frame our whole bounds, for debugging; note that we draw in only a portion of our bounds, and the rest is transparent
	//[[NSColor blackColor] set];
	//NSFrameRect(bounds);
	
#if 0
	// Debugging code: frame and fill our label rect without using NSBezierPath
	[[NSColor colorWithCalibratedHue:0.15 saturation:0.2 brightness:1.0 alpha:1.0] set];
	NSRectFill(labelRect);
	
	[[NSColor colorWithCalibratedHue:0.15 saturation:0.2 brightness:0.3 alpha:1.0] set];
	NSFrameRect(labelRect);
#else
	// Production code: use NSBezierPath to get a label that has a tag off of it
	NSBezierPath *bez = [NSBezierPath bezierPath];
	NSRect ilr = NSInsetRect(labelRect, 0.5, 0.5);	// inset by 0.5 to place us mid-pixel, so stroke looks good
	const double tagHeight = 5.0;
	
	if ([marker isLeftMarker])
	{
		// label rect with a diagonal tag down from the right edge
		[bez moveToPoint:NSMakePoint(ilr.origin.x, ilr.origin.y)];
		[bez relativeLineToPoint:NSMakePoint(ilr.size.width - tagHeight, 0)];
		[bez relativeLineToPoint:NSMakePoint(tagHeight, -tagHeight)];
		[bez relativeLineToPoint:NSMakePoint(0, ilr.size.height + tagHeight)];
		[bez relativeLineToPoint:NSMakePoint(-ilr.size.width, 0)];
		[bez closePath];
	}
	else
	{
		// label rect with a diagonal tag down from the right edge
		[bez moveToPoint:NSMakePoint(ilr.origin.x + ilr.size.width, ilr.origin.y)];
		[bez relativeLineToPoint:NSMakePoint(- (ilr.size.width - tagHeight), 0)];
		[bez relativeLineToPoint:NSMakePoint(-tagHeight, -tagHeight)];
		[bez relativeLineToPoint:NSMakePoint(0, ilr.size.height + tagHeight)];
		[bez relativeLineToPoint:NSMakePoint(ilr.size.width, 0)];
		[bez closePath];
	}
	
	[[NSColor colorWithCalibratedHue:0.15 saturation:0.2 brightness:1.0 alpha:1.0] set];
	[bez fill];
	
	[[NSColor blackColor] set];
	[bez setLineWidth:1.0];
	[bez stroke];
#endif
	
	[attrLabel drawAtPoint:NSMakePoint(labelRect.origin.x + 4, labelRect.origin.y + 1)];
	[attrLabel release];
}

- (BOOL)isOpaque
{
	return NO;
}

@end

@implementation SLiMSelectionMarker

// makes a new marker with no label and no tip point, not shown
+ (instancetype)new
{
	return [[[self class] alloc] initWithContentRect:NSMakeRect(0, 0, 150, 20) styleMask:NSBorderlessWindowMask backing:NSBackingStoreBuffered defer:YES];	// 150x20 should suffice, unless we change our font size...
}

- (instancetype)initWithContentRect:(NSRect)contentRect styleMask:(NSUInteger)aStyle backing:(NSBackingStoreType)bufferingType defer:(BOOL)flag
{
	if (self = [super initWithContentRect:contentRect styleMask:aStyle backing:bufferingType defer:flag])
	{
		[self setFloatingPanel:YES];
		[self setBecomesKeyOnlyIfNeeded:YES];
		[self setHasShadow:NO];
		[self setOpaque:NO];
		[self setBackgroundColor:[NSColor clearColor]];
		
		SLiMSelectionView *view = [[SLiMSelectionView alloc] initWithFrame:contentRect];
		
		[self setContentView:view];
		[view release];
		
		_tipPoint = NSMakePoint(round(contentRect.origin.x + contentRect.size.width / 2.0), contentRect.origin.y);
		_label = [@"10000000000" retain];
	}
	
	return self;
}

- (void)setLabel:(NSString *)label
{
	if (![_label isEqualToString:label])
	{
		[label retain];
		[_label release];
		_label = label;
		
		[[self contentView] setNeedsDisplay:YES];
	}
}

- (void)setTipPoint:(NSPoint)tipPoint
{
	if (!NSEqualPoints(_tipPoint, tipPoint))
	{
		NSPoint origin = [self frame].origin;
		
		origin.x += (tipPoint.x - _tipPoint.x);
		origin.y += (tipPoint.y - _tipPoint.y);
		
		_tipPoint = tipPoint;
		
		[self setFrameOrigin:origin];
	}
}

@end


@implementation NSScreen (SLiMWindowFrames)

+ (BOOL)visibleCandidateWindowFrame:(NSRect)candidateFrame
{
	NSArray *screens = [NSScreen screens];
	NSUInteger nScreens = [screens count];
	
	for (int i = 0; i < nScreens; ++i)
	{
		NSScreen *screen = [screens objectAtIndex:i];
		NSRect screenFrame = [screen visibleFrame];
		
		if (NSContainsRect(screenFrame, candidateFrame))
			return YES;
	}
	
	return NO;
}

@end


@implementation NSPopUpButton (SLiMSorting)

- (void)slimSortMenuItemsByTag
{
	NSMenu *menu = [self menu];
	int nItems = (int)[menu numberOfItems];
	
	// completely dumb bubble sort; not worth worrying about
	do
	{
		BOOL foundSwap = NO;
		
		for (int i = 0; i < nItems - 1; ++i)
		{
			NSMenuItem *firstItem = [menu itemAtIndex:i];
			NSMenuItem *secondItem = [menu itemAtIndex:i + 1];
			NSInteger firstTag = [firstItem tag];
			NSInteger secondTag = [secondItem tag];
			
			if (firstTag > secondTag)
			{
				[secondItem retain];
				[menu removeItemAtIndex:i + 1];
				[menu insertItem:secondItem atIndex:i];
				[secondItem release];
				
				foundSwap = YES;
			}
		}
		
		if (!foundSwap)
			break;
	}
	while (YES);
}

@end

@implementation NSPopUpButton (SLiMTinting)

- (void)slimSetTintColor:(NSColor *)tintColor
{
	[self setContentFilters:[NSArray array]];
	
	if (tintColor)
	{
		if (!self.layer)
			[self setWantsLayer:YES];
		
		CIFilter *tintFilter = [CIFilter filterWithName:@"CIColorMatrix"];
		
		if (tintFilter)
		{
			CGFloat redComponent = [tintColor redComponent];
			CGFloat greenComponent = [tintColor greenComponent];
			CGFloat blueComponent = [tintColor blueComponent];
			
			//NSLog(@"tintColor: redComponent == %f, greenComponent == %f, blueComponent == %f", redComponent, greenComponent, blueComponent);
			
			// The goal is to use CIColorMatrix to multiply color components so that white turns into tintColor; these vectors do that
			CIVector *rVector = [CIVector vectorWithX:redComponent Y:0.0 Z:0.0 W:0.0];
			CIVector *gVector = [CIVector vectorWithX:0.0 Y:greenComponent Z:0.0 W:0.0];
			CIVector *bVector = [CIVector vectorWithX:0.0 Y:0.0 Z:blueComponent W:0.0];
			
			[tintFilter setDefaults];
			[tintFilter setValue:rVector forKey:@"inputRVector"];
			[tintFilter setValue:gVector forKey:@"inputGVector"];
			[tintFilter setValue:bVector forKey:@"inputBVector"];
			
			[self setContentFilters:[NSArray arrayWithObject:tintFilter]];
		}
		else
		{
			NSLog(@"could not create [CIFilter filterWithName:@\"CIColorMatrix\"]");
		}
	}
	
	[self setNeedsDisplay];
	[self.layer setNeedsDisplay];
}

@end

#endif	// def SLIMGUI




































